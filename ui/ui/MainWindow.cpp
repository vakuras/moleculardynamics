///
/// Main window implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "stdafx.h"
#include "MainWindow.h"
#include "AboutWindow.h"

using namespace System;
using namespace System::Windows::Forms;
using namespace System::Drawing;
using namespace System::Diagnostics;

namespace ui 
{
	//retrieve and create an _instance of the form
	MainWindow ^ MainWindow::GetInstance()
	{
		if (_instance==nullptr)
		{
			_instance = gcnew MainWindow();
		}

		return _instance;
	}

	//constructor
	MainWindow::MainWindow(void)
	{
		_titleBar = gcnew ui::TitleBar();
		_titleBar->Visible = true;
		_titleBar->Location.X = 0;
		_titleBar->Location.Y = 0;
		_titleBar->Dock = System::Windows::Forms::DockStyle::Top;
		_titleBar->CloseClicked += gcnew EventHandler(this, &MainWindow::_titleBarEvent);
		_titleBar->SaveClicked += gcnew EventHandler(this, &MainWindow::_titleBarSave);
		_titleBar->SetTitle("Configuration");
		_titleBar->SaveSetVisible(true);
		_consoleBar = gcnew ui::TitleBar();
		_consoleBar->Visible = true;
		_consoleBar->Location.X = 0;
		_consoleBar->Location.Y = 0;
		_consoleBar->Dock = System::Windows::Forms::DockStyle::Top;
		_consoleBar->CloseClicked += gcnew EventHandler(this, &MainWindow::_consoleBarEvent);
		_consoleBar->SetTitle("Console");
		_compareBar = gcnew ui::TitleBar();
		_compareBar->Visible = true;
		_compareBar->Location.X = 0;
		_compareBar->Location.Y = 0;
		_compareBar->Dock = System::Windows::Forms::DockStyle::Top;
		_compareBar->CloseClicked += gcnew EventHandler(this, &MainWindow::_compareBarEvent);
		_compareBar->SetTitle("Compare");

		InitializeComponent();

		config = new configuration();

		//read configuration
		if (!ReadConfiguration("default.cfg", config))
		{
			UpdateConfiguration(config);
			WriteConfiguration("default.cfg", config);
		}

		//write to view
		UpdateView(config);

		_pnlCompare->Dock = System::Windows::Forms::DockStyle::Fill;
		_pnlPlot->Dock = System::Windows::Forms::DockStyle::Fill;
		_pnlResults->Dock = System::Windows::Forms::DockStyle::Fill;

		_plotWindow = gcnew PlotWindow(config);
		_plotWindow->Dock = System::Windows::Forms::DockStyle::Fill;
		_plotWindow->BackColor = _pnlCompareWindow->BackColor;
		_plotWindow->Parent = _pnlPlot;
		_plotWindow->BorderStyle = System::Windows::Forms::BorderStyle::None;
		_pnlPlot->Controls->Add(_plotWindow);

		_resultsWindow = gcnew ResultsWindow();
		_resultsWindow->Dock = System::Windows::Forms::DockStyle::Fill;
		_resultsWindow->BackColor = _pnlCompareWindow->BackColor;
		_resultsWindow->Parent = _pnlPlot;
		_resultsWindow->BorderStyle = System::Windows::Forms::BorderStyle::None;
		_pnlResults->Controls->Add(_resultsWindow);

		_pnlSplit->Panel1->Controls->Add(_titleBar);
		_split->Panel2->Controls->Add(_consoleBar);
		_pnlCompareWindow->Controls->Add(_compareBar);
	}

	//destructor
	MainWindow::~MainWindow()
	{
		//dispose of _instance
		_instance = nullptr;

		//update configuration structure with the data from the view
		UpdateConfiguration(config);

		//write it to the default configuration file
		WriteConfiguration("default.cfg", config);

		if (_components)
		{
			delete _components;
		}

		//if a _process is alive - kill it!
		if (_process!=nullptr && !_process->HasExited)
			_process->Kill();
	}

	System::Void MainWindow::_titleBarSave(System::Object^  , System::EventArgs^ )
	{
		//update configuration structure with the data from the view
		UpdateConfiguration(config);

		//write it to the default configuration file
		WriteConfiguration("default.cfg", config);
	}

	void MainWindow::UpdateConfiguration(configuration * config)
	{
		config->DT = float::Parse(_txtDT->Text);
		config->LennardJonesParticles = (int) _nudParticles->Value;
		config->ManyBodyParticles = (int) _nudMBParticles->Value;
		config->LennardJonesRCUT = float::Parse(_txtRCUT->Text);
		config->LennardJonesRS = float::Parse(_txtRS->Text);
		config->ManyBodyRCUT = float::Parse(_txtMBRCUT->Text);
		config->ManyBodyRS = float::Parse(_txtMBRS->Text);
		config->Temperature = float::Parse(_txtTemp->Text);
		config->Timesteps = (int) _nudTimesteps->Value;
		config->OutputTimesteps = (int) _nudOutputTimesteps->Value;
		config->useLennardJones = (int) _chkLJ->Checked;
		config->useManyBody = (int) _chkMB->Checked;
		config->lennardJonesBPP = (int) _chkLJBPP->Checked;
		config->manyBodyBPP = (int) _chkMBBPP->Checked;
		config->Fallback = (int) _chkFallback->Checked;
		config->Input = MakeCSTR(_txtInput->Text);
		config->Output = MakeCSTR(_txtOutput->Text);
		config->Filename = MakeCSTR(_txtResult->Text);
		config->UseCuda = _mnuCUDA->Checked;
		config->CudaBlocks = int::Parse(_tstCudaBlocks->Text);
		config->Debug = _chkDebug->Checked;
		config->vxcmFrom = float::Parse(_txtVxcmFrom->Text);
		config->vxcmTo = float::Parse(_txtVxcmTo->Text);
		config->vxcmStep = float::Parse(_txtVxcmStep->Text);
		config->AnimFile = MakeCSTR(_txtAnimFilename->Text);
		config->animts = (int) _nudAnimTS->Value;
	}

	void MainWindow::UpdateView(configuration * config)
	{
		_txtDT->Text = config->DT + "";
		_nudParticles->Value = config->LennardJonesParticles;
		_nudMBParticles->Value = config->ManyBodyParticles;
		_txtRCUT->Text = config->LennardJonesRCUT + "";
		_txtRS->Text = config->LennardJonesRS + "";
		_txtMBRCUT->Text = config->ManyBodyRCUT + "";
		_txtMBRS->Text = config->ManyBodyRS + "";
		_txtTemp->Text = config->Temperature + "";
		_nudTimesteps->Value = config->Timesteps;
		_nudOutputTimesteps->Value = config->OutputTimesteps;
		_nudAnimTS->Value = config->animts;

		_txtInput->Text = gcnew System::String(config->Input.c_str());
		_txtOutput->Text = gcnew System::String(config->Output.c_str());
		_txtResult->Text = gcnew System::String(config->Filename.c_str());
		_txtAnimFilename->Text = gcnew System::String(config->AnimFile.c_str());

		_mnuCUDA->Checked = config->UseCuda;
		_tsbCUDA->Checked = config->UseCuda;
		_chkLJ->Checked = config->useLennardJones;
		_chkMB->Checked = config->useManyBody;
		_chkLJBPP->Checked = config->lennardJonesBPP;
		_chkMBBPP->Checked = config->manyBodyBPP;
		_chkFallback->Checked = config->Fallback;
		_tstCudaBlocks->Text = config->CudaBlocks + "";
		_chkDebug->Checked = config->Debug;

		_txtVxcmFrom->Text = config->vxcmFrom + "";
		_txtVxcmTo->Text = config->vxcmTo + "";
		_txtVxcmStep->Text = config->vxcmStep + "";
	}

	System::Void MainWindow::_mnuLoadResults_Click(System::Object^  , System::EventArgs^  )
	{
		String^ filename = ChooseFile(false,"Open a results file", "Results files |*.dat");

		if (String::IsNullOrEmpty(filename))
			return;

		_mnuResult_Click(nullptr,nullptr);

		_resultsWindow->GetResults(filename);
	}

	System::Void MainWindow::_mnuLoadLog_Click(System::Object^  , System::EventArgs^  )
	{
		String^ filename = ChooseFile(false,"Open a results file", "Results files |*.txt");

		if (String::IsNullOrEmpty(filename))
			return;

		_split->Panel2Collapsed = false;

		_lstConsole->Items->Clear();

		cli::array<String^,1> ^ lines = System::IO::File::ReadAllLines(filename);

		for each(String ^ line in lines)
			_lstConsole->Items->Add(line);

		if (_lstConsole->Items->Count > 0)
			_lstConsole->SelectedIndex = _lstConsole->Items->Count - 1;
			
	}

	System::Void MainWindow::_mnuExit_Click(System::Object^  , System::EventArgs^  ) 
	{
		Application::Exit();
	}

	System::Void MainWindow::_mnuConsole_Click(System::Object^  , System::EventArgs^  ) 
	{
		_split->Panel2Collapsed = !_split->Panel2Collapsed;
	}

	System::Void MainWindow::_mnuConfiguration_Click(System::Object^  , System::EventArgs^  ) 
	{
		_pnlSplit->Panel1Collapsed = !_pnlSplit->Panel1Collapsed;
	}

	System::Void MainWindow::MainWindow_FormClosing(System::Object^  , System::Windows::Forms::FormClosingEventArgs^  ) 
	{
		Application::Exit();
	}

	System::Void MainWindow::_mnuCompare_Click(System::Object^  , System::EventArgs^  ) 
	{
		_pnlPlot->Visible = false;
		_plotWindow->Visible = false;
		_plotWindow->Destroy();
		_resultsWindow->Visible = false;
		_pnlCompare->Visible = true;
		_pnlResults->Visible = false;
	}

	System::Void MainWindow::_mnuCUDA_Click(System::Object^  , System::EventArgs^ ) 
	{
		if (_mnuCUDA->Checked)
		{
			_mnuCUDA->Checked = false;
			_tsbCUDA->Checked = false;
		}
		else
		{
			_mnuCUDA->Checked = true;
			_tsbCUDA->Checked = true;
		}
	}
	
	System::Void MainWindow::_tstCudaBlocks_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^ ) 
	{
		try
		{
			int::Parse(_tstCudaBlocks->Text);
		}
		catch(System::Exception ^)
		{
			_tstCudaBlocks->Text = config->CudaBlocks + ""; 
		}
	}

	System::Void MainWindow::_mnuExecute_Click(System::Object^  , System::EventArgs^  ) 
	{
		//update configuration file with the view data
		UpdateConfiguration(config);

		WriteConfiguration("default.cfg", config);

		//create a new _process
		_process = gcnew System::Diagnostics::Process();

		//depending on the cuda use checkbox - select the appropriate filename for the simulation
		if (_mnuCUDA->Checked)
			_process->StartInfo->FileName = Application::StartupPath + "\\cudamd.exe";
		else
			_process->StartInfo->FileName = Application::StartupPath + "\\cpumd.exe";

		_lstConsole->Items->Clear();
		
		_split->Panel2Collapsed = false;
		_pnlPlot->Visible = false;
		_plotWindow->Visible = false;
		_pnlCompare->Visible = false;
		_pnlResults->Visible = true;
		_plotWindow->Destroy();
		_resultsWindow->Visible = false;

		//_process configuration
		_process->StartInfo->Arguments = "\"" + Application::StartupPath + "\\default.cfg" + "\"";
		_process->StartInfo->UseShellExecute = false;
		_process->StartInfo->CreateNoWindow = true;
		_process->StartInfo->RedirectStandardError = true;
		_process->StartInfo->RedirectStandardOutput = true;
		_process->EnableRaisingEvents = true;
		_process->OutputDataReceived+= gcnew DataReceivedEventHandler(this, &ui::MainWindow::EventOutputHandler);
		_process->ErrorDataReceived+= gcnew DataReceivedEventHandler(this, &ui::MainWindow::EventOutputHandler);
		_process->Exited += gcnew EventHandler(this, &ui::MainWindow::EventProcessExited);

		//start the _process and check if the _process failed to start
		try
		{
			_process->Start();
		}
		catch(System::Exception ^ex)
		{
			MSGERR("Unable to run simulation!\n" + ex->Message);
			return;
		}

		//begin async redirection of standard output and standard error
		_process->BeginOutputReadLine();
		_process->BeginErrorReadLine();

		_pnlSplit->Panel1->Enabled = false;

		//update the status
		_tssStatus->Text = "Executing...";
	}

	void MainWindow::StatusReady()
	{
		//if it is not the current thread and the from requires to be invoked
		if (this->InvokeRequired)
		{
			//create an _instance of the delegate function and invoke it
			Invoke(gcnew MethodInvoker(this, &ui::MainWindow::StatusReady));
		}
		else //ready status
		{
			_tssStatus->Text = "Ready.";
		}
	}

	System::Void MainWindow::_mnuPlot_Click(System::Object^  , System::EventArgs^  ) 
	{
		_pnlPlot->Visible = true;
		_pnlCompare->Visible = false;
		_resultsWindow->Visible = false;
		_plotWindow->Visible = true;
	}

	System::Void MainWindow::EventOutputHandler(System::Object^ , System::Diagnostics::DataReceivedEventArgs^ outLine)
	{
		//if the data is null (empty)
		if (outLine->Data == nullptr)
			return;

		//otherwise - add the text to the console
		AddText(outLine->Data);
	}

	System::Void MainWindow::EventProcessExited(System::Object^ , System::EventArgs ^ )
	{
		StatusReady();
		PanelEnable();
	}

	void MainWindow::AddText(String ^ text)
	{
		//if it is not the current thread and _txtConsole requires to be invoked
		if (_lstConsole->InvokeRequired)
		{
			cli::array<Object ^>^ txtobj = gcnew cli::array<String ^>(1);
			txtobj[0] = text;
			//create an _instance of the delegate function and invoke it
			Invoke(gcnew AddTextCallback(this, &ui::MainWindow::AddText), txtobj);
		}
		else
		{
			_lstConsole->Items->Add(text);
			_lstConsole->SelectedIndex = _lstConsole->Items->Count - 1;
		}
	}

	//open or save dialog method - return filename on success or null on failure
	System::String ^ MainWindow::ChooseFile(bool save, System::String ^ title, System::String ^ filter)
	{
		if (save)
		{
			SaveFileDialog^ sfd = gcnew SaveFileDialog();
			sfd->InitialDirectory = Application::StartupPath;
			sfd->Filter = filter;
			sfd->FilterIndex = 1;
			sfd->RestoreDirectory = true;
			sfd->Title = title;

			if (sfd->ShowDialog() == System::Windows::Forms::DialogResult::OK)
				return sfd->FileName;
		}
		else
		{
			OpenFileDialog^ ofd = gcnew OpenFileDialog();
			ofd->InitialDirectory = Application::StartupPath;
			ofd->Filter = filter;
			ofd->FilterIndex = 1;
			ofd->RestoreDirectory = true;
			ofd->Title = title;

			if (ofd->ShowDialog() == System::Windows::Forms::DialogResult::OK)
				return ofd->FileName;
		}

		return nullptr;
	}

	void MainWindow::PanelEnable()
	{
		//if it is not the current thread and the from requires to be invoked
		if (_lstConsole->InvokeRequired)
		{
			//create an _instance of the delegate function and invoke it
			Invoke(gcnew MethodInvoker(this, &ui::MainWindow::PanelEnable));
		}
		else
		{
			_pnlSplit->Panel1->Enabled = true;
		}
	}

	System::Void MainWindow::_btnInput_Click(System::Object^  , System::EventArgs^  ) 
	{
		System::String ^ filename = ChooseFile(false, "Choose an input file", "All files |*.*");

		if (filename!=nullptr)
			_txtInput->Text = filename;
	}

	System::Void MainWindow::_btnOutput_Click(System::Object^  , System::EventArgs^  ) 
	{
		System::String ^ filename = ChooseFile(true, "Choose an output file", "All files |*.*");

		if (filename!=nullptr)
			_txtOutput->Text = filename;
	}

	System::Void MainWindow::_btnImport_Click(System::Object^  , System::EventArgs^  ) 
	{
		System::String ^ filename = ChooseFile(false,  "Import a configuration file", "Configuration file|*.cfg");

		if (filename!=nullptr)
		{
			string filenamestr;

			filenamestr = MakeCSTR(filename);

			ReadConfiguration(filenamestr, config);
			UpdateView(config);
		}
	}

	System::Void MainWindow::_btnExport_Click(System::Object^  , System::EventArgs^  ) 
	{
		System::String ^ filename = ChooseFile(true,  "Export a configuration file", "Configuration file|*.cfg");

		if (filename!=nullptr)
		{
			string filenamestr;

			filenamestr = MakeCSTR(filename);

			UpdateConfiguration(config);
			WriteConfiguration(filenamestr, config);
		}
	}

	//pop up to where menu
	System::Void MainWindow::_btnCompare1_Click(System::Object^  , System::EventArgs^  ) 
	{
		_contextButton = 1;

		_mnuContext->Show(_btnCompare1,0,0);
	}

	//pop up to where menu
	System::Void MainWindow::_btnCompare2_Click(System::Object^  , System::EventArgs^  )
	{
		_contextButton = 2;

		_mnuContext->Show(_btnCompare2,0,0);
	}

	System::Void MainWindow::_mnuCCompare_Click(System::Object^  , System::EventArgs^  ) 
	{
		if (_contextButton == 1)
			_txtInput1->Text = _txtInput->Text;
		else
			_txtInput2->Text = _txtOutput->Text;

		_mnuCompare_Click(this, nullptr);
	}


	System::Void MainWindow::_txtDT_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  ) 
	{
		try
		{
			float::Parse(_txtDT->Text);
		}
		catch(System::Exception ^)
		{
			_txtDT->Text = config->DT + ""; 
		}
	}

	System::Void MainWindow::_txtMBRS_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
	{
		try
		{
			float::Parse(_txtMBRS->Text);
		}
		catch(System::Exception ^)
		{
			_txtMBRS->Text = config->ManyBodyRS + ""; 
		}
	}

	System::Void MainWindow::_txtMBRCUT_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
	{
		try
		{
			float::Parse(_txtMBRCUT->Text);
		}
		catch(System::Exception ^)
		{
			_txtMBRCUT->Text = config->ManyBodyRCUT + ""; 
		}
	}

	System::Void MainWindow::_txtTemp_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  ) 
	{
		try
		{
			float::Parse(_txtTemp->Text);
		}
		catch(System::Exception ^)
		{
			_txtTemp->Text = config->Temperature + ""; 
		}
	}

	System::Void MainWindow::_txtEnergyLoss_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  ) 
	{
		try
		{
			float::Parse(_txtEnergyLoss->Text);
		}
		catch(System::Exception ^)
		{
			_txtEnergyLoss->Text = config->energyLoss + ""; 
		}
	}

	System::Void MainWindow::_txtVxcmFrom_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
	{
		try
		{
			float::Parse(_txtVxcmFrom->Text);
		}
		catch(System::Exception ^)
		{
			_txtVxcmFrom->Text = config->vxcmFrom + ""; 
		}
	}

	System::Void MainWindow::_txtVxcmTo_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
	{
		try
		{
			float::Parse(_txtVxcmTo->Text);
		}
		catch(System::Exception ^)
		{
			_txtVxcmTo->Text = config->vxcmTo + ""; 
		}
	}

	System::Void MainWindow::_txtVxcmStep_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
	{
		try
		{
			float::Parse(_txtVxcmStep->Text);
		}
		catch(System::Exception ^)
		{
			_txtVxcmStep->Text = config->vxcmStep + ""; 
		}
	}

	System::Void MainWindow::_txtRS_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  ) 
	{
		try
		{
			float::Parse(_txtRS->Text);
		}
		catch(System::Exception ^)
		{
			_txtRS->Text = config->LennardJonesRS + ""; 
		}
	}

	System::Void MainWindow::_txtRCUT_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  ) 
	{
		try
		{
			float::Parse(_txtRCUT->Text);
		}
		catch(System::Exception ^)
		{
			_txtRCUT->Text = config->LennardJonesRCUT + ""; 
		}
	}

	System::Void MainWindow::_mnuCPlot_Click(System::Object^  , System::EventArgs^  ) 
	{
		string str;

		if (_contextButton == 1)
			str = MakeCSTR(_txtInput->Text);
		else
			str = MakeCSTR(_txtOutput->Text);

		_plotWindow->LoadSimpleFile(str);

		_pnlPlot->Visible = true;
		_plotWindow->Visible = true;
		_pnlCompare->Visible = false;
		_pnlResults->Visible = false;
		_resultsWindow->Visible = false;
	}

	System::Void MainWindow::_pnlCompare_Resize(System::Object^  , System::EventArgs^  )
	{
		_pnlCompareWindow->Left = _pnlCompare->Width / 2 - _pnlCompareWindow->Width / 2;
		_pnlCompareWindow->Top = _pnlCompare->Height / 2 - _pnlCompareWindow->Height / 2;
	}

	System::Void MainWindow::_pnlResults_Resize(System::Object^  , System::EventArgs^  )
	{
		_resultsWindow->Left = _pnlResults->Width / 2 - _resultsWindow->Width / 2;
		_resultsWindow->Top = _pnlResults->Height / 2 - _resultsWindow->Height / 2;
	}

	System::Void MainWindow::_split_Panel2_Resize(System::Object^  , System::EventArgs^  )
	{
		_lstConsole->Left = 0;
		_lstConsole->Width = _split->Panel2->Width;
		_lstConsole->Top = _consoleBar->Bottom;
		_lstConsole->Height = _split->Panel2->Height - _consoleBar->Height;
	}

	System::Void MainWindow::_titleBarEvent(System::Object^ sender, System::EventArgs^ e)
	{
		_mnuConfiguration_Click(sender, e);
	}

	System::Void MainWindow::_consoleBarEvent(System::Object^ sender, System::EventArgs^ e)
	{
		_mnuConsole_Click(sender, e);
	}

	System::Void MainWindow::_compareBarEvent(System::Object^ sender, System::EventArgs^ e)
	{
		_pnlCompare->Visible = false;
	}

	System::Void MainWindow::_btnInput1_Click(System::Object^  , System::EventArgs^ )
	{
		System::String ^ filename = ChooseFile(false, "Choose the 1st input file", "Text files |*.txt");

		if (filename!=nullptr)
			_txtInput1->Text = filename;
	}
	
	System::Void MainWindow::_btnInput2_Click(System::Object^  , System::EventArgs^  )
	{
		System::String ^ filename = ChooseFile(false, "Choose the 2nd input file", "Text files |*.txt");

		if (filename!=nullptr)
			_txtInput2->Text = filename;
	}
	
	System::Void MainWindow::_btnCompare_Click(System::Object^  , System::EventArgs^  )
	{
		if (String::IsNullOrEmpty(_txtInput1->Text) || String::IsNullOrEmpty(_txtInput2->Text))
			return;

		UpdateConfiguration(config);

		int nop = (config->LennardJonesParticles + config->ManyBodyParticles);
		float * posArray1 = new float[nop * 4];
		float * posArray2 = new float[nop * 4];
		float * velocityArray1 = new float[nop * 3];
		float * velocityArray2 = new float[nop * 3];

		float px1;
		float px2;
		float pv1;
		float pv2;
		float dx;
		float dv;
		float hdx = 0;
		float hdv = 0;
		float pdx = 0;
		float pdv = 0;

		string input1;
		string input2;

		input1 = MakeCSTR(_txtInput1->Text);
		input2 = MakeCSTR(_txtInput2->Text);

		if (!ReadInput(input1.c_str(), nop, posArray1, velocityArray1) || !ReadInput(input2.c_str(), nop, posArray2, velocityArray2))
		{
			MSGERR("Unable to open input file!");
			return;
		}

		for(int i=0; i<nop; i++)
		{
			px1 = posArray1[i*4] * posArray1[i*4] + posArray1[i*4+1] * posArray1[i*4+1] + posArray1[i*4+2] * posArray1[i*4+2];
			px1 = (float) Math::Sqrt(px1);

			px2 = posArray2[i*4] * posArray2[i*4] + posArray2[i*4+1] * posArray2[i*4+1] + posArray2[i*4+2] * posArray2[i*4+2];
			px2 = (float) Math::Sqrt(px2);

			dx = px2 - px1;

			if (dx>hdx)
			{
				hdx = dx;
				pdx = (dx / px1) * 100;
			}

			pv1 = velocityArray1[i*3] * velocityArray1[i*3] + velocityArray1[i*3+1] * velocityArray1[i*3+1] + velocityArray1[i*3+2] * velocityArray1[i*3+2];
			pv1 = (float) Math::Sqrt(pv1);

			pv2 = velocityArray2[i*3] * velocityArray2[i*3] + velocityArray2[i*3+1] * velocityArray2[i*3+1] + velocityArray2[i*3+2] * velocityArray2[i*3+2];
			pv2 = (float) Math::Sqrt(pv2);

			dv = pv2 - pv1;

			if (dv>hdv)
			{
				hdv = dv;
				pdv = (dv / pv1) * 100;
			}
		}

		_lblCompare->Text = L"Largest δx: " + hdx + L" δx%: " + pdx + "%" +  L"\nLargest δv: " + hdv + L" δv%: " + pdv + "%";

		delete [] posArray1;
		delete [] posArray2;
		delete [] velocityArray1;
		delete [] velocityArray2;
	}
	
	System::Void MainWindow::_btnClipboard_Click(System::Object^  , System::EventArgs^  )
	{
		if (_lblCompare->Text==nullptr)
			return;

		Clipboard::SetText(_lblCompare->Text);
	}

	System::Void MainWindow::_mnuResult_Click(System::Object^  , System::EventArgs^  )
	{
		_pnlPlot->Visible = false;
		_plotWindow->Visible = false;
		_pnlCompare->Visible = false;
		_pnlResults->Visible = true;
		_plotWindow->Destroy();
		_resultsWindow->Visible = true;
	}

	///
	/// Open about window.
	///
	System::Void MainWindow::_mnuAbout_Click(System::Object^  , System::EventArgs^  )
	{
		(gcnew AboutWindow())->ShowDialog(this);
	}
}