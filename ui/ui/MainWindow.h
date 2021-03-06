///
/// Main window definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#pragma once

#ifndef MAINH
#define MAINH

using namespace System;
using namespace System::Windows::Forms;
using namespace System::Drawing;

#include "PlotWindow.h"
#include "ResultsWindow.h"
#include "TitleBar.h"

namespace ui 
{
	//Singleton main window form
	public ref class MainWindow : public System::Windows::Forms::Form
	{
	private:
		static MainWindow ^ _instance = nullptr; //holds the _instance of the form
		static System::Diagnostics::Process ^ _process;
		int _contextButton;
		configuration * config; //configuration struct
		ui::PlotWindow ^ _plotWindow;
		ui::ResultsWindow ^ _resultsWindow;
		System::ComponentModel::IContainer^  _components;
		System::Windows::Forms::MenuStrip^  _mnuMainMenu;
		System::Windows::Forms::StatusStrip^  _stsStatusbar;
		System::Windows::Forms::ToolStripMenuItem^  _mnuFile;
		System::Windows::Forms::ToolStripMenuItem^  _mnuExit;
		System::Windows::Forms::ToolStripMenuItem^  _mnuWindow;
		System::Windows::Forms::ToolStripMenuItem^  _mnuConfiguration;
		System::Windows::Forms::ToolStripMenuItem^  _mnuConsole;
		System::Windows::Forms::ToolStripStatusLabel^  _tssStatus;
		System::Windows::Forms::ToolStripMenuItem^  _mnuCompare;
		System::Windows::Forms::ToolStrip^  _tsToolbar;
		System::Windows::Forms::ToolStripButton^  _tsbCUDA;
		System::Windows::Forms::ToolStripButton^  _tsbExecute;
		System::Windows::Forms::ToolStripButton^  _tsbCompare;
		System::Windows::Forms::ToolStripSeparator^  _tsSeperator;
		System::Windows::Forms::ToolStripButton^  _tsbConfiguration;
		System::Windows::Forms::ToolStripButton^  _tsbConsole;
		System::Windows::Forms::ToolStripLabel^  _tslCudaBlocks;
		System::Windows::Forms::ToolStripTextBox^  _tstCudaBlocks;
		System::Windows::Forms::ToolStripMenuItem^  _mnuCUDA;
		System::Windows::Forms::ToolStripMenuItem^  _mnuExecute;
		System::Windows::Forms::ToolStripSeparator^  _mnuSeperator1;
		System::Windows::Forms::ToolStripMenuItem^  _mnuPlot;
		System::Windows::Forms::ToolStripButton^  _tsbPlot;
		System::Windows::Forms::SplitContainer^  _split;
		System::Windows::Forms::SplitContainer^  _pnlSplit;
		System::Windows::Forms::TextBox^  _txtDT;
		System::Windows::Forms::Label^  _lblDT;
		System::Windows::Forms::Label^  _lblTimesteps;
		System::Windows::Forms::NumericUpDown^  _nudTimesteps;
		System::Windows::Forms::NumericUpDown^  _nudParticles;
		System::Windows::Forms::Label^  _lblParticles;
		System::Windows::Forms::TextBox^  _txtEnergyLoss;
		System::Windows::Forms::TextBox^  _txtTemp;
		System::Windows::Forms::Label^  _lblEnergyLoss;
		System::Windows::Forms::Label^  _lblTemp;
		System::Windows::Forms::TextBox^  _txtRCUT;
		System::Windows::Forms::TextBox^  _txtRS;
		System::Windows::Forms::Label^  _lblRCUT;
		System::Windows::Forms::Label^  _lblRS;
		System::Windows::Forms::Button^  _btnCompare1;
		System::Windows::Forms::Button^  _btnCompare2;
		System::Windows::Forms::Button^  _btnOutput;
		System::Windows::Forms::Button^  _btnInput;
		System::Windows::Forms::Label^  _lblOutput;
		System::Windows::Forms::TextBox^  _txtInput;
		System::Windows::Forms::Label^  _lblInput;
		System::Windows::Forms::TextBox^  _txtOutput;
		System::Windows::Forms::Label^  _lblResult;
		System::Windows::Forms::TextBox^  _txtResult;
		System::Windows::Forms::Button^  _btnImport;
		System::Windows::Forms::Button^  _btnExport;
		System::Windows::Forms::ContextMenuStrip^  _mnuContext;
		System::Windows::Forms::ToolStripMenuItem^  _mnuCCompare;
		System::Windows::Forms::ToolStripMenuItem^  _mnuCPlot;
		System::Windows::Forms::Panel^  _pnlPlot;
		System::Windows::Forms::Panel^  _pnlResults;
		System::Windows::Forms::Panel^  _pnlCompare;
		System::Windows::Forms::Panel^  _pnlCompareWindow;
		System::Windows::Forms::Button^  _btnInput2;
		System::Windows::Forms::Button^  _btnInput1;
		System::Windows::Forms::TextBox^  _txtInput1;
		System::Windows::Forms::Label^  _lblInput1;
		System::Windows::Forms::TextBox^  _txtInput2;
		System::Windows::Forms::Label^  _lblInput2;
		System::Windows::Forms::Label^  _lblCompare;
		System::Windows::Forms::Button^  _btnClipboard;
		System::Windows::Forms::Button^  _btnCompare;
		System::Windows::Forms::ToolStripButton^  _tsbResult;
		System::Windows::Forms::ToolStripMenuItem^  _mnuResult;
		System::Windows::Forms::ToolStripMenuItem^  _mnuHelp;
		System::Windows::Forms::ToolStripMenuItem^  _mnuAbout;
		System::ComponentModel::IContainer^  components;
		TitleBar^ _titleBar;
		TitleBar^ _consoleBar;
		System::Windows::Forms::TextBox^  _txtMBRCUT;
		System::Windows::Forms::TextBox^  _txtMBRS;
		System::Windows::Forms::Label^  _lblMRCUT;
		System::Windows::Forms::Label^  _lblMBRS;
		System::Windows::Forms::Label^  _lblMBParticles;
		System::Windows::Forms::NumericUpDown^  _nudMBParticles;
		System::Windows::Forms::CheckBox^  _chkDebug;
		System::Windows::Forms::Label^  _lblOutputTimeSteps;
		System::Windows::Forms::NumericUpDown^  _nudOutputTimesteps;
		System::Windows::Forms::CheckBox^  _chkMB;
		System::Windows::Forms::CheckBox^  _chkLJ;
		System::Windows::Forms::Panel^  _pnlConfigPanel;
		System::Windows::Forms::CheckBox^  _chkFallback;
		System::Windows::Forms::CheckBox^  _chkMBBPP;
		System::Windows::Forms::CheckBox^  _chkLJBPP;
		System::Windows::Forms::TextBox^  _txtVxcmStep;
		System::Windows::Forms::TextBox^  _txtVxcmTo;
		System::Windows::Forms::Label^  label1;
		System::Windows::Forms::Label^  label2;
		System::Windows::Forms::TextBox^  _txtVxcmFrom;
		System::Windows::Forms::Label^  label3;
		System::Windows::Forms::ListBox^  _lstConsole;
		System::Windows::Forms::Label^  _lblAnimTS;
		System::Windows::Forms::NumericUpDown^  _nudAnimTS;
		System::Windows::Forms::Label^  _lblAnimFilename;
		System::Windows::Forms::TextBox^  _txtAnimFilename;
		System::Windows::Forms::ToolStripButton^  _tsLoadLog;
		System::Windows::Forms::ToolStripButton^  _tsLoadResults;
		System::Windows::Forms::ToolStripSeparator^  _tsSeperator2;
		System::Windows::Forms::ToolStripMenuItem^  _mnuLoadLog;
		System::Windows::Forms::ToolStripMenuItem^  _mnuLoadResults;
		System::Windows::Forms::ToolStripSeparator^  _mnuSeperator2;
		TitleBar^ _compareBar;

	private:
		System::Void _mnuLoadResults_Click(System::Object^  , System::EventArgs^  );
		System::Void _mnuLoadLog_Click(System::Object^  , System::EventArgs^  );
		System::Void _txtVxcmFrom_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e);
		System::Void _txtVxcmTo_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e);
		System::Void _txtVxcmStep_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e);
		System::Void _mnuResult_Click(System::Object^  , System::EventArgs^  );
		System::Void EventOutputHandler(System::Object^ , System::Diagnostics::DataReceivedEventArgs^ outLine);
		System::Void EventProcessExited(System::Object^ , System::EventArgs ^ );
		System::Void _mnuExit_Click(System::Object^  , System::EventArgs^  );
		System::Void _mnuConsole_Click(System::Object^  , System::EventArgs^  );
		System::Void _mnuConfiguration_Click(System::Object^  , System::EventArgs^  );
		System::Void MainWindow_FormClosing(System::Object^  , System::Windows::Forms::FormClosingEventArgs^  );
		System::Void _mnuCompare_Click(System::Object^  , System::EventArgs^  ); 
		System::Void _mnuCUDA_Click(System::Object^  , System::EventArgs^ ); 
		System::Void _tstCudaBlocks_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^ );
		System::Void _mnuExecute_Click(System::Object^  , System::EventArgs^  );
		System::Void _mnuCPlot_Click(System::Object^  , System::EventArgs^  );
		System::Void _btnInput_Click(System::Object^  , System::EventArgs^  ); 
		System::Void _btnOutput_Click(System::Object^  , System::EventArgs^  );
		System::Void _btnImport_Click(System::Object^  , System::EventArgs^  ); 
		System::Void _btnExport_Click(System::Object^  , System::EventArgs^  ); 
		System::Void _btnCompare1_Click(System::Object^  , System::EventArgs^  ); 
		System::Void _btnCompare2_Click(System::Object^  , System::EventArgs^  );
		System::Void _mnuCCompare_Click(System::Object^  , System::EventArgs^  );
		System::Void _txtDT_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  );
		System::Void _txtTemp_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  );
		System::Void _txtEnergyLoss_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  );
		System::Void _txtRS_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  );
		System::Void _txtRCUT_Validating(System::Object^  , System::ComponentModel::CancelEventArgs^  );
		System::Void _mnuPlot_Click(System::Object^  , System::EventArgs^  );
		System::Void _split_Panel2_Resize(System::Object^  , System::EventArgs^  );
		System::Void _btnInput1_Click(System::Object^  , System::EventArgs^ );
		System::Void _btnInput2_Click(System::Object^  , System::EventArgs^  );
		System::Void _btnCompare_Click(System::Object^  , System::EventArgs^  );
		System::Void _btnClipboard_Click(System::Object^  , System::EventArgs^  );
		System::Void _pnlCompare_Resize(System::Object^  , System::EventArgs^  );
		System::Void _pnlResults_Resize(System::Object^  , System::EventArgs^  );
		System::Void _mnuAbout_Click(System::Object^  , System::EventArgs^  );
		System::Void _titleBarEvent(System::Object^  , System::EventArgs^ );
		System::Void _titleBarSave(System::Object^  , System::EventArgs^ );
		System::Void _consoleBarEvent(System::Object^  , System::EventArgs^ );
		System::Void _compareBarEvent(System::Object^  , System::EventArgs^ );
		System::Void _txtMBRS_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e);
		System::Void _txtMBRCUT_Validating(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e);
		System::String ^ ChooseFile(bool save, System::String ^ , System::String ^ );
		void PanelEnable();
		void StatusReady();
		void UpdateConfiguration(configuration * );
		void UpdateView(configuration * );
		void AddText(String ^ );
		delegate void AddTextCallback(String ^);
	public:
		static MainWindow ^ GetInstance();

	protected:
		MainWindow(void);
		~MainWindow();

#pragma region Windows Form Designer generated code
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(MainWindow::typeid));
			this->_mnuMainMenu = (gcnew System::Windows::Forms::MenuStrip());
			this->_mnuFile = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuCUDA = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuExecute = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuSeperator1 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->_mnuLoadLog = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuLoadResults = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuSeperator2 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->_mnuExit = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuWindow = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuConfiguration = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuCompare = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuConsole = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuPlot = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuResult = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuHelp = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuAbout = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_stsStatusbar = (gcnew System::Windows::Forms::StatusStrip());
			this->_tssStatus = (gcnew System::Windows::Forms::ToolStripStatusLabel());
			this->_tsToolbar = (gcnew System::Windows::Forms::ToolStrip());
			this->_tsbConfiguration = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsbConsole = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsbCompare = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsbPlot = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsbResult = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsSeperator = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->_tsLoadLog = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsLoadResults = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tsSeperator2 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->_tsbCUDA = (gcnew System::Windows::Forms::ToolStripButton());
			this->_tslCudaBlocks = (gcnew System::Windows::Forms::ToolStripLabel());
			this->_tstCudaBlocks = (gcnew System::Windows::Forms::ToolStripTextBox());
			this->_tsbExecute = (gcnew System::Windows::Forms::ToolStripButton());
			this->_split = (gcnew System::Windows::Forms::SplitContainer());
			this->_pnlSplit = (gcnew System::Windows::Forms::SplitContainer());
			this->_pnlConfigPanel = (gcnew System::Windows::Forms::Panel());
			this->_lblAnimTS = (gcnew System::Windows::Forms::Label());
			this->_nudAnimTS = (gcnew System::Windows::Forms::NumericUpDown());
			this->_lblAnimFilename = (gcnew System::Windows::Forms::Label());
			this->_txtAnimFilename = (gcnew System::Windows::Forms::TextBox());
			this->_txtVxcmStep = (gcnew System::Windows::Forms::TextBox());
			this->_txtVxcmTo = (gcnew System::Windows::Forms::TextBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->_txtVxcmFrom = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->_chkFallback = (gcnew System::Windows::Forms::CheckBox());
			this->_chkMBBPP = (gcnew System::Windows::Forms::CheckBox());
			this->_chkLJBPP = (gcnew System::Windows::Forms::CheckBox());
			this->_chkMB = (gcnew System::Windows::Forms::CheckBox());
			this->_chkLJ = (gcnew System::Windows::Forms::CheckBox());
			this->_lblOutputTimeSteps = (gcnew System::Windows::Forms::Label());
			this->_nudOutputTimesteps = (gcnew System::Windows::Forms::NumericUpDown());
			this->_chkDebug = (gcnew System::Windows::Forms::CheckBox());
			this->_lblMBParticles = (gcnew System::Windows::Forms::Label());
			this->_nudMBParticles = (gcnew System::Windows::Forms::NumericUpDown());
			this->_txtMBRCUT = (gcnew System::Windows::Forms::TextBox());
			this->_txtMBRS = (gcnew System::Windows::Forms::TextBox());
			this->_lblMRCUT = (gcnew System::Windows::Forms::Label());
			this->_lblMBRS = (gcnew System::Windows::Forms::Label());
			this->_btnExport = (gcnew System::Windows::Forms::Button());
			this->_btnImport = (gcnew System::Windows::Forms::Button());
			this->_lblResult = (gcnew System::Windows::Forms::Label());
			this->_txtResult = (gcnew System::Windows::Forms::TextBox());
			this->_btnCompare1 = (gcnew System::Windows::Forms::Button());
			this->_btnCompare2 = (gcnew System::Windows::Forms::Button());
			this->_btnOutput = (gcnew System::Windows::Forms::Button());
			this->_btnInput = (gcnew System::Windows::Forms::Button());
			this->_lblOutput = (gcnew System::Windows::Forms::Label());
			this->_txtInput = (gcnew System::Windows::Forms::TextBox());
			this->_lblInput = (gcnew System::Windows::Forms::Label());
			this->_txtOutput = (gcnew System::Windows::Forms::TextBox());
			this->_txtRCUT = (gcnew System::Windows::Forms::TextBox());
			this->_txtRS = (gcnew System::Windows::Forms::TextBox());
			this->_lblRCUT = (gcnew System::Windows::Forms::Label());
			this->_lblRS = (gcnew System::Windows::Forms::Label());
			this->_txtEnergyLoss = (gcnew System::Windows::Forms::TextBox());
			this->_txtTemp = (gcnew System::Windows::Forms::TextBox());
			this->_lblEnergyLoss = (gcnew System::Windows::Forms::Label());
			this->_lblTemp = (gcnew System::Windows::Forms::Label());
			this->_txtDT = (gcnew System::Windows::Forms::TextBox());
			this->_lblDT = (gcnew System::Windows::Forms::Label());
			this->_lblTimesteps = (gcnew System::Windows::Forms::Label());
			this->_nudTimesteps = (gcnew System::Windows::Forms::NumericUpDown());
			this->_nudParticles = (gcnew System::Windows::Forms::NumericUpDown());
			this->_lblParticles = (gcnew System::Windows::Forms::Label());
			this->_pnlPlot = (gcnew System::Windows::Forms::Panel());
			this->_pnlResults = (gcnew System::Windows::Forms::Panel());
			this->_pnlCompare = (gcnew System::Windows::Forms::Panel());
			this->_pnlCompareWindow = (gcnew System::Windows::Forms::Panel());
			this->_lblCompare = (gcnew System::Windows::Forms::Label());
			this->_btnClipboard = (gcnew System::Windows::Forms::Button());
			this->_btnCompare = (gcnew System::Windows::Forms::Button());
			this->_lblInput2 = (gcnew System::Windows::Forms::Label());
			this->_btnInput2 = (gcnew System::Windows::Forms::Button());
			this->_btnInput1 = (gcnew System::Windows::Forms::Button());
			this->_txtInput1 = (gcnew System::Windows::Forms::TextBox());
			this->_lblInput1 = (gcnew System::Windows::Forms::Label());
			this->_txtInput2 = (gcnew System::Windows::Forms::TextBox());
			this->_lstConsole = (gcnew System::Windows::Forms::ListBox());
			this->_mnuContext = (gcnew System::Windows::Forms::ContextMenuStrip(this->components));
			this->_mnuCCompare = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuCPlot = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->_mnuMainMenu->SuspendLayout();
			this->_stsStatusbar->SuspendLayout();
			this->_tsToolbar->SuspendLayout();
			this->_split->Panel1->SuspendLayout();
			this->_split->Panel2->SuspendLayout();
			this->_split->SuspendLayout();
			this->_pnlSplit->Panel1->SuspendLayout();
			this->_pnlSplit->Panel2->SuspendLayout();
			this->_pnlSplit->SuspendLayout();
			this->_pnlConfigPanel->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudAnimTS))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudOutputTimesteps))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudMBParticles))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudTimesteps))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudParticles))->BeginInit();
			this->_pnlCompare->SuspendLayout();
			this->_pnlCompareWindow->SuspendLayout();
			this->_mnuContext->SuspendLayout();
			this->SuspendLayout();
			// 
			// _mnuMainMenu
			// 
			this->_mnuMainMenu->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(3) {this->_mnuFile, this->_mnuWindow, 
				this->_mnuHelp});
			this->_mnuMainMenu->Location = System::Drawing::Point(0, 0);
			this->_mnuMainMenu->Name = L"_mnuMainMenu";
			this->_mnuMainMenu->Size = System::Drawing::Size(840, 24);
			this->_mnuMainMenu->TabIndex = 1;
			// 
			// _mnuFile
			// 
			this->_mnuFile->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(7) {this->_mnuCUDA, this->_mnuExecute, 
				this->_mnuSeperator1, this->_mnuLoadLog, this->_mnuLoadResults, this->_mnuSeperator2, this->_mnuExit});
			this->_mnuFile->Name = L"_mnuFile";
			this->_mnuFile->Size = System::Drawing::Size(37, 20);
			this->_mnuFile->Text = L"&File";
			// 
			// _mnuCUDA
			// 
			this->_mnuCUDA->Checked = true;
			this->_mnuCUDA->CheckState = System::Windows::Forms::CheckState::Checked;
			this->_mnuCUDA->Name = L"_mnuCUDA";
			this->_mnuCUDA->Size = System::Drawing::Size(152, 22);
			this->_mnuCUDA->Text = L"&Use CUDA";
			this->_mnuCUDA->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCUDA_Click);
			// 
			// _mnuExecute
			// 
			this->_mnuExecute->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuExecute.Image")));
			this->_mnuExecute->Name = L"_mnuExecute";
			this->_mnuExecute->Size = System::Drawing::Size(152, 22);
			this->_mnuExecute->Text = L"&Execute";
			this->_mnuExecute->Click += gcnew System::EventHandler(this, &MainWindow::_mnuExecute_Click);
			// 
			// _mnuSeperator1
			// 
			this->_mnuSeperator1->Name = L"_mnuSeperator1";
			this->_mnuSeperator1->Size = System::Drawing::Size(149, 6);
			// 
			// _mnuLoadLog
			// 
			this->_mnuLoadLog->Name = L"_mnuLoadLog";
			this->_mnuLoadLog->Size = System::Drawing::Size(152, 22);
			this->_mnuLoadLog->Text = L"Load &Log";
			this->_mnuLoadLog->Click += gcnew System::EventHandler(this, &MainWindow::_mnuLoadLog_Click);
			// 
			// _mnuLoadResults
			// 
			this->_mnuLoadResults->Name = L"_mnuLoadResults";
			this->_mnuLoadResults->Size = System::Drawing::Size(152, 22);
			this->_mnuLoadResults->Text = L"Load &Results";
			this->_mnuLoadResults->Click += gcnew System::EventHandler(this, &MainWindow::_mnuLoadResults_Click);
			// 
			// _mnuSeperator2
			// 
			this->_mnuSeperator2->Name = L"_mnuSeperator2";
			this->_mnuSeperator2->Size = System::Drawing::Size(149, 6);
			// 
			// _mnuExit
			// 
			this->_mnuExit->Name = L"_mnuExit";
			this->_mnuExit->Size = System::Drawing::Size(152, 22);
			this->_mnuExit->Text = L"E&xit";
			this->_mnuExit->Click += gcnew System::EventHandler(this, &MainWindow::_mnuExit_Click);
			// 
			// _mnuWindow
			// 
			this->_mnuWindow->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(5) {this->_mnuConfiguration, 
				this->_mnuCompare, this->_mnuConsole, this->_mnuPlot, this->_mnuResult});
			this->_mnuWindow->Name = L"_mnuWindow";
			this->_mnuWindow->Size = System::Drawing::Size(63, 20);
			this->_mnuWindow->Text = L"&Window";
			// 
			// _mnuConfiguration
			// 
			this->_mnuConfiguration->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuConfiguration.Image")));
			this->_mnuConfiguration->Name = L"_mnuConfiguration";
			this->_mnuConfiguration->Size = System::Drawing::Size(162, 22);
			this->_mnuConfiguration->Text = L"&Configuration";
			this->_mnuConfiguration->Click += gcnew System::EventHandler(this, &MainWindow::_mnuConfiguration_Click);
			// 
			// _mnuCompare
			// 
			this->_mnuCompare->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuCompare.Image")));
			this->_mnuCompare->Name = L"_mnuCompare";
			this->_mnuCompare->Size = System::Drawing::Size(162, 22);
			this->_mnuCompare->Text = L"C&ompare (δx δv)";
			this->_mnuCompare->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCompare_Click);
			// 
			// _mnuConsole
			// 
			this->_mnuConsole->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuConsole.Image")));
			this->_mnuConsole->Name = L"_mnuConsole";
			this->_mnuConsole->Size = System::Drawing::Size(162, 22);
			this->_mnuConsole->Text = L"Con&sole";
			this->_mnuConsole->Click += gcnew System::EventHandler(this, &MainWindow::_mnuConsole_Click);
			// 
			// _mnuPlot
			// 
			this->_mnuPlot->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuPlot.Image")));
			this->_mnuPlot->Name = L"_mnuPlot";
			this->_mnuPlot->Size = System::Drawing::Size(162, 22);
			this->_mnuPlot->Text = L"&Plot";
			this->_mnuPlot->Click += gcnew System::EventHandler(this, &MainWindow::_mnuPlot_Click);
			// 
			// _mnuResult
			// 
			this->_mnuResult->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_mnuResult.Image")));
			this->_mnuResult->Name = L"_mnuResult";
			this->_mnuResult->Size = System::Drawing::Size(162, 22);
			this->_mnuResult->Text = L"Result Graph";
			this->_mnuResult->Click += gcnew System::EventHandler(this, &MainWindow::_mnuResult_Click);
			// 
			// _mnuHelp
			// 
			this->_mnuHelp->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {this->_mnuAbout});
			this->_mnuHelp->Name = L"_mnuHelp";
			this->_mnuHelp->Size = System::Drawing::Size(44, 20);
			this->_mnuHelp->Text = L"&Help";
			// 
			// _mnuAbout
			// 
			this->_mnuAbout->Name = L"_mnuAbout";
			this->_mnuAbout->Size = System::Drawing::Size(107, 22);
			this->_mnuAbout->Text = L"&About";
			this->_mnuAbout->Click += gcnew System::EventHandler(this, &MainWindow::_mnuAbout_Click);
			// 
			// _stsStatusbar
			// 
			this->_stsStatusbar->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {this->_tssStatus});
			this->_stsStatusbar->Location = System::Drawing::Point(0, 749);
			this->_stsStatusbar->Name = L"_stsStatusbar";
			this->_stsStatusbar->Size = System::Drawing::Size(840, 22);
			this->_stsStatusbar->TabIndex = 2;
			// 
			// _tssStatus
			// 
			this->_tssStatus->Name = L"_tssStatus";
			this->_tssStatus->Size = System::Drawing::Size(42, 17);
			this->_tssStatus->Text = L"Ready.";
			this->_tssStatus->TextDirection = System::Windows::Forms::ToolStripTextDirection::Horizontal;
			// 
			// _tsToolbar
			// 
			this->_tsToolbar->GripStyle = System::Windows::Forms::ToolStripGripStyle::Hidden;
			this->_tsToolbar->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(13) {this->_tsbConfiguration, 
				this->_tsbConsole, this->_tsbCompare, this->_tsbPlot, this->_tsbResult, this->_tsSeperator, this->_tsLoadLog, this->_tsLoadResults, 
				this->_tsSeperator2, this->_tsbCUDA, this->_tslCudaBlocks, this->_tstCudaBlocks, this->_tsbExecute});
			this->_tsToolbar->Location = System::Drawing::Point(0, 24);
			this->_tsToolbar->Name = L"_tsToolbar";
			this->_tsToolbar->Size = System::Drawing::Size(840, 25);
			this->_tsToolbar->TabIndex = 4;
			// 
			// _tsbConfiguration
			// 
			this->_tsbConfiguration->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbConfiguration.Image")));
			this->_tsbConfiguration->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbConfiguration->Name = L"_tsbConfiguration";
			this->_tsbConfiguration->Size = System::Drawing::Size(101, 22);
			this->_tsbConfiguration->Text = L"Configuration";
			this->_tsbConfiguration->Click += gcnew System::EventHandler(this, &MainWindow::_mnuConfiguration_Click);
			// 
			// _tsbConsole
			// 
			this->_tsbConsole->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbConsole.Image")));
			this->_tsbConsole->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbConsole->Name = L"_tsbConsole";
			this->_tsbConsole->Size = System::Drawing::Size(70, 22);
			this->_tsbConsole->Text = L"Console";
			this->_tsbConsole->Click += gcnew System::EventHandler(this, &MainWindow::_mnuConsole_Click);
			// 
			// _tsbCompare
			// 
			this->_tsbCompare->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbCompare.Image")));
			this->_tsbCompare->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbCompare->Name = L"_tsbCompare";
			this->_tsbCompare->Size = System::Drawing::Size(76, 22);
			this->_tsbCompare->Text = L"Compare";
			this->_tsbCompare->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCompare_Click);
			// 
			// _tsbPlot
			// 
			this->_tsbPlot->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbPlot.Image")));
			this->_tsbPlot->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbPlot->Name = L"_tsbPlot";
			this->_tsbPlot->Size = System::Drawing::Size(48, 22);
			this->_tsbPlot->Text = L"Plot";
			this->_tsbPlot->Click += gcnew System::EventHandler(this, &MainWindow::_mnuPlot_Click);
			// 
			// _tsbResult
			// 
			this->_tsbResult->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbResult.Image")));
			this->_tsbResult->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbResult->Name = L"_tsbResult";
			this->_tsbResult->Size = System::Drawing::Size(94, 22);
			this->_tsbResult->Text = L"Result Graph";
			this->_tsbResult->Click += gcnew System::EventHandler(this, &MainWindow::_mnuResult_Click);
			// 
			// _tsSeperator
			// 
			this->_tsSeperator->Name = L"_tsSeperator";
			this->_tsSeperator->Size = System::Drawing::Size(6, 25);
			// 
			// _tsLoadLog
			// 
			this->_tsLoadLog->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsLoadLog.Image")));
			this->_tsLoadLog->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsLoadLog->Name = L"_tsLoadLog";
			this->_tsLoadLog->Size = System::Drawing::Size(76, 22);
			this->_tsLoadLog->Text = L"Load Log";
			this->_tsLoadLog->Click += gcnew System::EventHandler(this, &MainWindow::_mnuLoadLog_Click);
			// 
			// _tsLoadResults
			// 
			this->_tsLoadResults->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsLoadResults.Image")));
			this->_tsLoadResults->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsLoadResults->Name = L"_tsLoadResults";
			this->_tsLoadResults->Size = System::Drawing::Size(93, 22);
			this->_tsLoadResults->Text = L"Load Results";
			this->_tsLoadResults->Click += gcnew System::EventHandler(this, &MainWindow::_mnuLoadResults_Click);
			// 
			// _tsSeperator2
			// 
			this->_tsSeperator2->Name = L"_tsSeperator2";
			this->_tsSeperator2->Size = System::Drawing::Size(6, 25);
			// 
			// _tsbCUDA
			// 
			this->_tsbCUDA->Checked = true;
			this->_tsbCUDA->CheckOnClick = true;
			this->_tsbCUDA->CheckState = System::Windows::Forms::CheckState::Checked;
			this->_tsbCUDA->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbCUDA.Image")));
			this->_tsbCUDA->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbCUDA->Name = L"_tsbCUDA";
			this->_tsbCUDA->Size = System::Drawing::Size(81, 22);
			this->_tsbCUDA->Text = L"Use CUDA";
			this->_tsbCUDA->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCUDA_Click);
			// 
			// _tslCudaBlocks
			// 
			this->_tslCudaBlocks->Name = L"_tslCudaBlocks";
			this->_tslCudaBlocks->Size = System::Drawing::Size(79, 22);
			this->_tslCudaBlocks->Text = L"CUDA Blocks:";
			// 
			// _tstCudaBlocks
			// 
			this->_tstCudaBlocks->MaxLength = 6;
			this->_tstCudaBlocks->Name = L"_tstCudaBlocks";
			this->_tstCudaBlocks->Size = System::Drawing::Size(40, 25);
			this->_tstCudaBlocks->Text = L"8";
			this->_tstCudaBlocks->TextBoxTextAlign = System::Windows::Forms::HorizontalAlignment::Right;
			this->_tstCudaBlocks->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_tstCudaBlocks_Validating);
			// 
			// _tsbExecute
			// 
			this->_tsbExecute->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"_tsbExecute.Image")));
			this->_tsbExecute->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->_tsbExecute->Name = L"_tsbExecute";
			this->_tsbExecute->Size = System::Drawing::Size(67, 20);
			this->_tsbExecute->Text = L"Execute";
			this->_tsbExecute->Click += gcnew System::EventHandler(this, &MainWindow::_mnuExecute_Click);
			// 
			// _split
			// 
			this->_split->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
			this->_split->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_split->Location = System::Drawing::Point(0, 49);
			this->_split->Name = L"_split";
			this->_split->Orientation = System::Windows::Forms::Orientation::Horizontal;
			// 
			// _split.Panel1
			// 
			this->_split->Panel1->Controls->Add(this->_pnlSplit);
			this->_split->Panel1MinSize = 300;
			// 
			// _split.Panel2
			// 
			this->_split->Panel2->Controls->Add(this->_lstConsole);
			this->_split->Panel2->Resize += gcnew System::EventHandler(this, &MainWindow::_split_Panel2_Resize);
			this->_split->Panel2Collapsed = true;
			this->_split->Size = System::Drawing::Size(840, 700);
			this->_split->SplitterDistance = 500;
			this->_split->TabIndex = 5;
			// 
			// _pnlSplit
			// 
			this->_pnlSplit->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
			this->_pnlSplit->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_pnlSplit->FixedPanel = System::Windows::Forms::FixedPanel::Panel1;
			this->_pnlSplit->IsSplitterFixed = true;
			this->_pnlSplit->Location = System::Drawing::Point(0, 0);
			this->_pnlSplit->Name = L"_pnlSplit";
			// 
			// _pnlSplit.Panel1
			// 
			this->_pnlSplit->Panel1->Controls->Add(this->_pnlConfigPanel);
			this->_pnlSplit->Panel1MinSize = 200;
			// 
			// _pnlSplit.Panel2
			// 
			this->_pnlSplit->Panel2->BackColor = System::Drawing::SystemColors::AppWorkspace;
			this->_pnlSplit->Panel2->Controls->Add(this->_pnlPlot);
			this->_pnlSplit->Panel2->Controls->Add(this->_pnlResults);
			this->_pnlSplit->Panel2->Controls->Add(this->_pnlCompare);
			this->_pnlSplit->Size = System::Drawing::Size(840, 700);
			this->_pnlSplit->SplitterDistance = 215;
			this->_pnlSplit->TabIndex = 0;
			// 
			// _pnlConfigPanel
			// 
			this->_pnlConfigPanel->AutoScroll = true;
			this->_pnlConfigPanel->Controls->Add(this->_lblAnimTS);
			this->_pnlConfigPanel->Controls->Add(this->_nudAnimTS);
			this->_pnlConfigPanel->Controls->Add(this->_lblAnimFilename);
			this->_pnlConfigPanel->Controls->Add(this->_txtAnimFilename);
			this->_pnlConfigPanel->Controls->Add(this->_txtVxcmStep);
			this->_pnlConfigPanel->Controls->Add(this->_txtVxcmTo);
			this->_pnlConfigPanel->Controls->Add(this->label1);
			this->_pnlConfigPanel->Controls->Add(this->label2);
			this->_pnlConfigPanel->Controls->Add(this->_txtVxcmFrom);
			this->_pnlConfigPanel->Controls->Add(this->label3);
			this->_pnlConfigPanel->Controls->Add(this->_chkFallback);
			this->_pnlConfigPanel->Controls->Add(this->_chkMBBPP);
			this->_pnlConfigPanel->Controls->Add(this->_chkLJBPP);
			this->_pnlConfigPanel->Controls->Add(this->_chkMB);
			this->_pnlConfigPanel->Controls->Add(this->_chkLJ);
			this->_pnlConfigPanel->Controls->Add(this->_lblOutputTimeSteps);
			this->_pnlConfigPanel->Controls->Add(this->_nudOutputTimesteps);
			this->_pnlConfigPanel->Controls->Add(this->_chkDebug);
			this->_pnlConfigPanel->Controls->Add(this->_lblMBParticles);
			this->_pnlConfigPanel->Controls->Add(this->_nudMBParticles);
			this->_pnlConfigPanel->Controls->Add(this->_txtMBRCUT);
			this->_pnlConfigPanel->Controls->Add(this->_txtMBRS);
			this->_pnlConfigPanel->Controls->Add(this->_lblMRCUT);
			this->_pnlConfigPanel->Controls->Add(this->_lblMBRS);
			this->_pnlConfigPanel->Controls->Add(this->_btnExport);
			this->_pnlConfigPanel->Controls->Add(this->_btnImport);
			this->_pnlConfigPanel->Controls->Add(this->_lblResult);
			this->_pnlConfigPanel->Controls->Add(this->_txtResult);
			this->_pnlConfigPanel->Controls->Add(this->_btnCompare1);
			this->_pnlConfigPanel->Controls->Add(this->_btnCompare2);
			this->_pnlConfigPanel->Controls->Add(this->_btnOutput);
			this->_pnlConfigPanel->Controls->Add(this->_btnInput);
			this->_pnlConfigPanel->Controls->Add(this->_lblOutput);
			this->_pnlConfigPanel->Controls->Add(this->_txtInput);
			this->_pnlConfigPanel->Controls->Add(this->_lblInput);
			this->_pnlConfigPanel->Controls->Add(this->_txtOutput);
			this->_pnlConfigPanel->Controls->Add(this->_txtRCUT);
			this->_pnlConfigPanel->Controls->Add(this->_txtRS);
			this->_pnlConfigPanel->Controls->Add(this->_lblRCUT);
			this->_pnlConfigPanel->Controls->Add(this->_lblRS);
			this->_pnlConfigPanel->Controls->Add(this->_txtEnergyLoss);
			this->_pnlConfigPanel->Controls->Add(this->_txtTemp);
			this->_pnlConfigPanel->Controls->Add(this->_lblEnergyLoss);
			this->_pnlConfigPanel->Controls->Add(this->_lblTemp);
			this->_pnlConfigPanel->Controls->Add(this->_txtDT);
			this->_pnlConfigPanel->Controls->Add(this->_lblDT);
			this->_pnlConfigPanel->Controls->Add(this->_lblTimesteps);
			this->_pnlConfigPanel->Controls->Add(this->_nudTimesteps);
			this->_pnlConfigPanel->Controls->Add(this->_nudParticles);
			this->_pnlConfigPanel->Controls->Add(this->_lblParticles);
			this->_pnlConfigPanel->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_pnlConfigPanel->Location = System::Drawing::Point(0, 0);
			this->_pnlConfigPanel->Name = L"_pnlConfigPanel";
			this->_pnlConfigPanel->Size = System::Drawing::Size(213, 698);
			this->_pnlConfigPanel->TabIndex = 69;
			// 
			// _lblAnimTS
			// 
			this->_lblAnimTS->AutoSize = true;
			this->_lblAnimTS->Location = System::Drawing::Point(5, 117);
			this->_lblAnimTS->Name = L"_lblAnimTS";
			this->_lblAnimTS->Size = System::Drawing::Size(107, 13);
			this->_lblAnimTS->TabIndex = 81;
			this->_lblAnimTS->Text = L"Animation Timesteps:";
			// 
			// _nudAnimTS
			// 
			this->_nudAnimTS->Location = System::Drawing::Point(116, 115);
			this->_nudAnimTS->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2147483647, 0, 0, 0});
			this->_nudAnimTS->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) {1, 0, 0, System::Int32::MinValue});
			this->_nudAnimTS->Name = L"_nudAnimTS";
			this->_nudAnimTS->Size = System::Drawing::Size(69, 20);
			this->_nudAnimTS->TabIndex = 80;
			this->_nudAnimTS->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {100, 0, 0, 0});
			// 
			// _lblAnimFilename
			// 
			this->_lblAnimFilename->AutoSize = true;
			this->_lblAnimFilename->Location = System::Drawing::Point(9, 729);
			this->_lblAnimFilename->Name = L"_lblAnimFilename";
			this->_lblAnimFilename->Size = System::Drawing::Size(122, 13);
			this->_lblAnimFilename->TabIndex = 78;
			this->_lblAnimFilename->Text = L"Animation data filename:";
			// 
			// _txtAnimFilename
			// 
			this->_txtAnimFilename->Location = System::Drawing::Point(12, 749);
			this->_txtAnimFilename->Name = L"_txtAnimFilename";
			this->_txtAnimFilename->Size = System::Drawing::Size(173, 20);
			this->_txtAnimFilename->TabIndex = 79;
			this->_txtAnimFilename->Text = L"anim";
			// 
			// _txtVxcmStep
			// 
			this->_txtVxcmStep->Location = System::Drawing::Point(116, 375);
			this->_txtVxcmStep->Name = L"_txtVxcmStep";
			this->_txtVxcmStep->Size = System::Drawing::Size(70, 20);
			this->_txtVxcmStep->TabIndex = 77;
			this->_txtVxcmStep->Text = L"0.01";
			this->_txtVxcmStep->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtVxcmStep_Validating);
			// 
			// _txtVxcmTo
			// 
			this->_txtVxcmTo->Location = System::Drawing::Point(117, 349);
			this->_txtVxcmTo->Name = L"_txtVxcmTo";
			this->_txtVxcmTo->Size = System::Drawing::Size(69, 20);
			this->_txtVxcmTo->TabIndex = 76;
			this->_txtVxcmTo->Text = L"0.12";
			this->_txtVxcmTo->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtVxcmTo_Validating);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(45, 378);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(65, 13);
			this->label1->TabIndex = 75;
			this->label1->Text = L"VXCM Step:";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(55, 352);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(56, 13);
			this->label2->TabIndex = 74;
			this->label2->Text = L"VXCM To:";
			// 
			// _txtVxcmFrom
			// 
			this->_txtVxcmFrom->Location = System::Drawing::Point(116, 323);
			this->_txtVxcmFrom->Name = L"_txtVxcmFrom";
			this->_txtVxcmFrom->Size = System::Drawing::Size(70, 20);
			this->_txtVxcmFrom->TabIndex = 73;
			this->_txtVxcmFrom->Text = L"0.05";
			this->_txtVxcmFrom->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtVxcmFrom_Validating);
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(45, 326);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(66, 13);
			this->label3->TabIndex = 72;
			this->label3->Text = L"VXCM From:";
			// 
			// _chkFallback
			// 
			this->_chkFallback->AutoSize = true;
			this->_chkFallback->Location = System::Drawing::Point(120, 521);
			this->_chkFallback->Name = L"_chkFallback";
			this->_chkFallback->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkFallback->Size = System::Drawing::Size(66, 17);
			this->_chkFallback->TabIndex = 71;
			this->_chkFallback->Text = L"Fallback";
			this->_chkFallback->UseVisualStyleBackColor = true;
			// 
			// _chkMBBPP
			// 
			this->_chkMBBPP->AutoSize = true;
			this->_chkMBBPP->Location = System::Drawing::Point(86, 498);
			this->_chkMBBPP->Name = L"_chkMBBPP";
			this->_chkMBBPP->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkMBBPP->Size = System::Drawing::Size(100, 17);
			this->_chkMBBPP->TabIndex = 70;
			this->_chkMBBPP->Text = L"ManyBody BPP";
			this->_chkMBBPP->UseVisualStyleBackColor = true;
			// 
			// _chkLJBPP
			// 
			this->_chkLJBPP->AutoSize = true;
			this->_chkLJBPP->Location = System::Drawing::Point(66, 475);
			this->_chkLJBPP->Name = L"_chkLJBPP";
			this->_chkLJBPP->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkLJBPP->Size = System::Drawing::Size(120, 17);
			this->_chkLJBPP->TabIndex = 69;
			this->_chkLJBPP->Text = L"Lennard Jones BPP";
			this->_chkLJBPP->UseVisualStyleBackColor = true;
			// 
			// _chkMB
			// 
			this->_chkMB->AutoSize = true;
			this->_chkMB->Location = System::Drawing::Point(88, 452);
			this->_chkMB->Name = L"_chkMB";
			this->_chkMB->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkMB->Size = System::Drawing::Size(98, 17);
			this->_chkMB->TabIndex = 68;
			this->_chkMB->Text = L"Use ManyBody";
			this->_chkMB->UseVisualStyleBackColor = true;
			// 
			// _chkLJ
			// 
			this->_chkLJ->AutoSize = true;
			this->_chkLJ->Location = System::Drawing::Point(68, 429);
			this->_chkLJ->Name = L"_chkLJ";
			this->_chkLJ->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkLJ->Size = System::Drawing::Size(118, 17);
			this->_chkLJ->TabIndex = 67;
			this->_chkLJ->Text = L"Use Lennard Jones";
			this->_chkLJ->UseVisualStyleBackColor = true;
			// 
			// _lblOutputTimeSteps
			// 
			this->_lblOutputTimeSteps->AutoSize = true;
			this->_lblOutputTimeSteps->Location = System::Drawing::Point(16, 91);
			this->_lblOutputTimeSteps->Name = L"_lblOutputTimeSteps";
			this->_lblOutputTimeSteps->Size = System::Drawing::Size(93, 13);
			this->_lblOutputTimeSteps->TabIndex = 66;
			this->_lblOutputTimeSteps->Text = L"Output Timesteps:";
			// 
			// _nudOutputTimesteps
			// 
			this->_nudOutputTimesteps->Location = System::Drawing::Point(116, 89);
			this->_nudOutputTimesteps->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2147483647, 0, 0, 0});
			this->_nudOutputTimesteps->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) {1, 0, 0, System::Int32::MinValue});
			this->_nudOutputTimesteps->Name = L"_nudOutputTimesteps";
			this->_nudOutputTimesteps->Size = System::Drawing::Size(69, 20);
			this->_nudOutputTimesteps->TabIndex = 65;
			this->_nudOutputTimesteps->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {100, 0, 0, 0});
			// 
			// _chkDebug
			// 
			this->_chkDebug->AutoSize = true;
			this->_chkDebug->Location = System::Drawing::Point(128, 406);
			this->_chkDebug->Name = L"_chkDebug";
			this->_chkDebug->RightToLeft = System::Windows::Forms::RightToLeft::Yes;
			this->_chkDebug->Size = System::Drawing::Size(58, 17);
			this->_chkDebug->TabIndex = 64;
			this->_chkDebug->Text = L"Debug";
			this->_chkDebug->UseVisualStyleBackColor = true;
			// 
			// _lblMBParticles
			// 
			this->_lblMBParticles->AutoSize = true;
			this->_lblMBParticles->Location = System::Drawing::Point(7, 39);
			this->_lblMBParticles->Name = L"_lblMBParticles";
			this->_lblMBParticles->Size = System::Drawing::Size(102, 13);
			this->_lblMBParticles->TabIndex = 63;
			this->_lblMBParticles->Text = L"ManyBody particles:";
			// 
			// _nudMBParticles
			// 
			this->_nudMBParticles->Location = System::Drawing::Point(116, 37);
			this->_nudMBParticles->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2147483647, 0, 0, 0});
			this->_nudMBParticles->Name = L"_nudMBParticles";
			this->_nudMBParticles->Size = System::Drawing::Size(69, 20);
			this->_nudMBParticles->TabIndex = 62;
			this->_nudMBParticles->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {28, 0, 0, 0});
			// 
			// _txtMBRCUT
			// 
			this->_txtMBRCUT->Location = System::Drawing::Point(116, 297);
			this->_txtMBRCUT->Name = L"_txtMBRCUT";
			this->_txtMBRCUT->Size = System::Drawing::Size(70, 20);
			this->_txtMBRCUT->TabIndex = 61;
			this->_txtMBRCUT->Text = L"20";
			this->_txtMBRCUT->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtMBRCUT_Validating);
			// 
			// _txtMBRS
			// 
			this->_txtMBRS->Location = System::Drawing::Point(117, 271);
			this->_txtMBRS->Name = L"_txtMBRS";
			this->_txtMBRS->Size = System::Drawing::Size(69, 20);
			this->_txtMBRS->TabIndex = 60;
			this->_txtMBRS->Text = L"22";
			this->_txtMBRS->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtMBRS_Validating);
			// 
			// _lblMRCUT
			// 
			this->_lblMRCUT->AutoSize = true;
			this->_lblMRCUT->Location = System::Drawing::Point(50, 300);
			this->_lblMRCUT->Name = L"_lblMRCUT";
			this->_lblMRCUT->Size = System::Drawing::Size(59, 13);
			this->_lblMRCUT->TabIndex = 59;
			this->_lblMRCUT->Text = L"MB RCUT:";
			// 
			// _lblMBRS
			// 
			this->_lblMBRS->AutoSize = true;
			this->_lblMBRS->Location = System::Drawing::Point(65, 274);
			this->_lblMBRS->Name = L"_lblMBRS";
			this->_lblMBRS->Size = System::Drawing::Size(44, 13);
			this->_lblMBRS->TabIndex = 58;
			this->_lblMBRS->Text = L"MB RS:";
			// 
			// _btnExport
			// 
			this->_btnExport->Location = System::Drawing::Point(12, 838);
			this->_btnExport->Name = L"_btnExport";
			this->_btnExport->Size = System::Drawing::Size(173, 28);
			this->_btnExport->TabIndex = 57;
			this->_btnExport->Text = L"Export Configuration";
			this->_btnExport->UseVisualStyleBackColor = true;
			this->_btnExport->Click += gcnew System::EventHandler(this, &MainWindow::_btnExport_Click);
			// 
			// _btnImport
			// 
			this->_btnImport->Location = System::Drawing::Point(12, 804);
			this->_btnImport->Name = L"_btnImport";
			this->_btnImport->Size = System::Drawing::Size(174, 28);
			this->_btnImport->TabIndex = 56;
			this->_btnImport->Text = L"Import Configuration";
			this->_btnImport->UseVisualStyleBackColor = true;
			this->_btnImport->Click += gcnew System::EventHandler(this, &MainWindow::_btnImport_Click);
			// 
			// _lblResult
			// 
			this->_lblResult->AutoSize = true;
			this->_lblResult->Location = System::Drawing::Point(9, 678);
			this->_lblResult->Name = L"_lblResult";
			this->_lblResult->Size = System::Drawing::Size(75, 13);
			this->_lblResult->TabIndex = 48;
			this->_lblResult->Text = L"Data filename:";
			// 
			// _txtResult
			// 
			this->_txtResult->Location = System::Drawing::Point(12, 694);
			this->_txtResult->Name = L"_txtResult";
			this->_txtResult->Size = System::Drawing::Size(173, 20);
			this->_txtResult->TabIndex = 49;
			this->_txtResult->Text = L"result";
			// 
			// _btnCompare1
			// 
			this->_btnCompare1->Location = System::Drawing::Point(158, 559);
			this->_btnCompare1->Name = L"_btnCompare1";
			this->_btnCompare1->Size = System::Drawing::Size(28, 23);
			this->_btnCompare1->TabIndex = 47;
			this->_btnCompare1->Text = L"»";
			this->_btnCompare1->UseVisualStyleBackColor = true;
			this->_btnCompare1->Click += gcnew System::EventHandler(this, &MainWindow::_btnCompare1_Click);
			// 
			// _btnCompare2
			// 
			this->_btnCompare2->Location = System::Drawing::Point(158, 610);
			this->_btnCompare2->Name = L"_btnCompare2";
			this->_btnCompare2->Size = System::Drawing::Size(28, 23);
			this->_btnCompare2->TabIndex = 46;
			this->_btnCompare2->Text = L"»";
			this->_btnCompare2->UseVisualStyleBackColor = true;
			this->_btnCompare2->Click += gcnew System::EventHandler(this, &MainWindow::_btnCompare2_Click);
			// 
			// _btnOutput
			// 
			this->_btnOutput->Location = System::Drawing::Point(124, 610);
			this->_btnOutput->Name = L"_btnOutput";
			this->_btnOutput->Size = System::Drawing::Size(28, 23);
			this->_btnOutput->TabIndex = 45;
			this->_btnOutput->Text = L"..";
			this->_btnOutput->UseVisualStyleBackColor = true;
			this->_btnOutput->Click += gcnew System::EventHandler(this, &MainWindow::_btnOutput_Click);
			// 
			// _btnInput
			// 
			this->_btnInput->Location = System::Drawing::Point(124, 559);
			this->_btnInput->Name = L"_btnInput";
			this->_btnInput->Size = System::Drawing::Size(28, 23);
			this->_btnInput->TabIndex = 44;
			this->_btnInput->Text = L"..";
			this->_btnInput->UseVisualStyleBackColor = true;
			this->_btnInput->Click += gcnew System::EventHandler(this, &MainWindow::_btnInput_Click);
			// 
			// _lblOutput
			// 
			this->_lblOutput->AutoSize = true;
			this->_lblOutput->Location = System::Drawing::Point(9, 619);
			this->_lblOutput->Name = L"_lblOutput";
			this->_lblOutput->Size = System::Drawing::Size(58, 13);
			this->_lblOutput->TabIndex = 43;
			this->_lblOutput->Text = L"Output file:";
			// 
			// _txtInput
			// 
			this->_txtInput->Location = System::Drawing::Point(12, 584);
			this->_txtInput->Name = L"_txtInput";
			this->_txtInput->Size = System::Drawing::Size(174, 20);
			this->_txtInput->TabIndex = 42;
			this->_txtInput->Text = L"input";
			// 
			// _lblInput
			// 
			this->_lblInput->AutoSize = true;
			this->_lblInput->Location = System::Drawing::Point(9, 564);
			this->_lblInput->Name = L"_lblInput";
			this->_lblInput->Size = System::Drawing::Size(50, 13);
			this->_lblInput->TabIndex = 41;
			this->_lblInput->Text = L"Input file:";
			// 
			// _txtOutput
			// 
			this->_txtOutput->Location = System::Drawing::Point(12, 639);
			this->_txtOutput->Name = L"_txtOutput";
			this->_txtOutput->Size = System::Drawing::Size(174, 20);
			this->_txtOutput->TabIndex = 40;
			this->_txtOutput->Text = L"output";
			// 
			// _txtRCUT
			// 
			this->_txtRCUT->Location = System::Drawing::Point(116, 245);
			this->_txtRCUT->Name = L"_txtRCUT";
			this->_txtRCUT->Size = System::Drawing::Size(70, 20);
			this->_txtRCUT->TabIndex = 36;
			this->_txtRCUT->Text = L"20";
			this->_txtRCUT->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtRCUT_Validating);
			// 
			// _txtRS
			// 
			this->_txtRS->Location = System::Drawing::Point(117, 219);
			this->_txtRS->Name = L"_txtRS";
			this->_txtRS->Size = System::Drawing::Size(69, 20);
			this->_txtRS->TabIndex = 35;
			this->_txtRS->Text = L"22";
			this->_txtRS->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtRS_Validating);
			// 
			// _lblRCUT
			// 
			this->_lblRCUT->AutoSize = true;
			this->_lblRCUT->Location = System::Drawing::Point(55, 248);
			this->_lblRCUT->Name = L"_lblRCUT";
			this->_lblRCUT->Size = System::Drawing::Size(54, 13);
			this->_lblRCUT->TabIndex = 34;
			this->_lblRCUT->Text = L"LJ RCUT:";
			// 
			// _lblRS
			// 
			this->_lblRS->AutoSize = true;
			this->_lblRS->Location = System::Drawing::Point(70, 222);
			this->_lblRS->Name = L"_lblRS";
			this->_lblRS->Size = System::Drawing::Size(39, 13);
			this->_lblRS->TabIndex = 33;
			this->_lblRS->Text = L"LJ RS:";
			// 
			// _txtEnergyLoss
			// 
			this->_txtEnergyLoss->Location = System::Drawing::Point(116, 193);
			this->_txtEnergyLoss->Name = L"_txtEnergyLoss";
			this->_txtEnergyLoss->Size = System::Drawing::Size(69, 20);
			this->_txtEnergyLoss->TabIndex = 30;
			this->_txtEnergyLoss->Text = L"0.2";
			this->_txtEnergyLoss->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtEnergyLoss_Validating);
			// 
			// _txtTemp
			// 
			this->_txtTemp->Location = System::Drawing::Point(117, 167);
			this->_txtTemp->Name = L"_txtTemp";
			this->_txtTemp->Size = System::Drawing::Size(69, 20);
			this->_txtTemp->TabIndex = 29;
			this->_txtTemp->Text = L"216.0";
			this->_txtTemp->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtTemp_Validating);
			// 
			// _lblEnergyLoss
			// 
			this->_lblEnergyLoss->AutoSize = true;
			this->_lblEnergyLoss->Location = System::Drawing::Point(43, 196);
			this->_lblEnergyLoss->Name = L"_lblEnergyLoss";
			this->_lblEnergyLoss->Size = System::Drawing::Size(68, 13);
			this->_lblEnergyLoss->TabIndex = 28;
			this->_lblEnergyLoss->Text = L"Energy Loss:";
			// 
			// _lblTemp
			// 
			this->_lblTemp->AutoSize = true;
			this->_lblTemp->Location = System::Drawing::Point(39, 170);
			this->_lblTemp->Name = L"_lblTemp";
			this->_lblTemp->Size = System::Drawing::Size(70, 13);
			this->_lblTemp->TabIndex = 27;
			this->_lblTemp->Text = L"Temperature:";
			// 
			// _txtDT
			// 
			this->_txtDT->Location = System::Drawing::Point(117, 141);
			this->_txtDT->Name = L"_txtDT";
			this->_txtDT->Size = System::Drawing::Size(69, 20);
			this->_txtDT->TabIndex = 14;
			this->_txtDT->Text = L"0.5";
			this->_txtDT->Validating += gcnew System::ComponentModel::CancelEventHandler(this, &MainWindow::_txtDT_Validating);
			// 
			// _lblDT
			// 
			this->_lblDT->AutoSize = true;
			this->_lblDT->Location = System::Drawing::Point(90, 144);
			this->_lblDT->Name = L"_lblDT";
			this->_lblDT->Size = System::Drawing::Size(19, 13);
			this->_lblDT->TabIndex = 13;
			this->_lblDT->Text = L"δt:";
			// 
			// _lblTimesteps
			// 
			this->_lblTimesteps->AutoSize = true;
			this->_lblTimesteps->Location = System::Drawing::Point(51, 65);
			this->_lblTimesteps->Name = L"_lblTimesteps";
			this->_lblTimesteps->Size = System::Drawing::Size(58, 13);
			this->_lblTimesteps->TabIndex = 12;
			this->_lblTimesteps->Text = L"Timesteps:";
			// 
			// _nudTimesteps
			// 
			this->_nudTimesteps->Location = System::Drawing::Point(116, 63);
			this->_nudTimesteps->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2147483647, 0, 0, 0});
			this->_nudTimesteps->Name = L"_nudTimesteps";
			this->_nudTimesteps->Size = System::Drawing::Size(69, 20);
			this->_nudTimesteps->TabIndex = 11;
			this->_nudTimesteps->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {10000, 0, 0, 0});
			// 
			// _nudParticles
			// 
			this->_nudParticles->Location = System::Drawing::Point(116, 11);
			this->_nudParticles->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2147483647, 0, 0, 0});
			this->_nudParticles->Name = L"_nudParticles";
			this->_nudParticles->Size = System::Drawing::Size(69, 20);
			this->_nudParticles->TabIndex = 10;
			this->_nudParticles->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {125, 0, 0, 0});
			// 
			// _lblParticles
			// 
			this->_lblParticles->AutoSize = true;
			this->_lblParticles->Location = System::Drawing::Point(41, 13);
			this->_lblParticles->Name = L"_lblParticles";
			this->_lblParticles->Size = System::Drawing::Size(66, 13);
			this->_lblParticles->TabIndex = 9;
			this->_lblParticles->Text = L"L-J particles:";
			// 
			// _pnlPlot
			// 
			this->_pnlPlot->BackColor = System::Drawing::SystemColors::AppWorkspace;
			this->_pnlPlot->Location = System::Drawing::Point(378, 526);
			this->_pnlPlot->Name = L"_pnlPlot";
			this->_pnlPlot->Size = System::Drawing::Size(115, 107);
			this->_pnlPlot->TabIndex = 3;
			this->_pnlPlot->Visible = false;
			// 
			// _pnlResults
			// 
			this->_pnlResults->BackColor = System::Drawing::SystemColors::AppWorkspace;
			this->_pnlResults->Location = System::Drawing::Point(57, 133);
			this->_pnlResults->Name = L"_pnlResults";
			this->_pnlResults->Size = System::Drawing::Size(517, 107);
			this->_pnlResults->TabIndex = 2;
			this->_pnlResults->Visible = false;
			this->_pnlResults->Resize += gcnew System::EventHandler(this, &MainWindow::_pnlResults_Resize);
			// 
			// _pnlCompare
			// 
			this->_pnlCompare->BackColor = System::Drawing::SystemColors::AppWorkspace;
			this->_pnlCompare->Controls->Add(this->_pnlCompareWindow);
			this->_pnlCompare->Location = System::Drawing::Point(57, 246);
			this->_pnlCompare->Name = L"_pnlCompare";
			this->_pnlCompare->Size = System::Drawing::Size(517, 245);
			this->_pnlCompare->TabIndex = 1;
			this->_pnlCompare->Visible = false;
			this->_pnlCompare->Resize += gcnew System::EventHandler(this, &MainWindow::_pnlCompare_Resize);
			// 
			// _pnlCompareWindow
			// 
			this->_pnlCompareWindow->BackColor = System::Drawing::SystemColors::Control;
			this->_pnlCompareWindow->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
			this->_pnlCompareWindow->Controls->Add(this->_lblCompare);
			this->_pnlCompareWindow->Controls->Add(this->_btnClipboard);
			this->_pnlCompareWindow->Controls->Add(this->_btnCompare);
			this->_pnlCompareWindow->Controls->Add(this->_lblInput2);
			this->_pnlCompareWindow->Controls->Add(this->_btnInput2);
			this->_pnlCompareWindow->Controls->Add(this->_btnInput1);
			this->_pnlCompareWindow->Controls->Add(this->_txtInput1);
			this->_pnlCompareWindow->Controls->Add(this->_lblInput1);
			this->_pnlCompareWindow->Controls->Add(this->_txtInput2);
			this->_pnlCompareWindow->Location = System::Drawing::Point(24, 22);
			this->_pnlCompareWindow->Name = L"_pnlCompareWindow";
			this->_pnlCompareWindow->Size = System::Drawing::Size(306, 175);
			this->_pnlCompareWindow->TabIndex = 1;
			// 
			// _lblCompare
			// 
			this->_lblCompare->Location = System::Drawing::Point(14, 118);
			this->_lblCompare->Name = L"_lblCompare";
			this->_lblCompare->Size = System::Drawing::Size(273, 42);
			this->_lblCompare->TabIndex = 57;
			this->_lblCompare->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// _btnClipboard
			// 
			this->_btnClipboard->Location = System::Drawing::Point(121, 85);
			this->_btnClipboard->Name = L"_btnClipboard";
			this->_btnClipboard->Size = System::Drawing::Size(166, 23);
			this->_btnClipboard->TabIndex = 56;
			this->_btnClipboard->Text = L"Copy To Clipboard";
			this->_btnClipboard->UseVisualStyleBackColor = true;
			this->_btnClipboard->Click += gcnew System::EventHandler(this, &MainWindow::_btnClipboard_Click);
			// 
			// _btnCompare
			// 
			this->_btnCompare->Location = System::Drawing::Point(17, 85);
			this->_btnCompare->Name = L"_btnCompare";
			this->_btnCompare->Size = System::Drawing::Size(98, 23);
			this->_btnCompare->TabIndex = 55;
			this->_btnCompare->Text = L"Compare";
			this->_btnCompare->UseVisualStyleBackColor = true;
			this->_btnCompare->Click += gcnew System::EventHandler(this, &MainWindow::_btnCompare_Click);
			// 
			// _lblInput2
			// 
			this->_lblInput2->AutoSize = true;
			this->_lblInput2->Location = System::Drawing::Point(14, 62);
			this->_lblInput2->Name = L"_lblInput2";
			this->_lblInput2->Size = System::Drawing::Size(59, 13);
			this->_lblInput2->TabIndex = 54;
			this->_lblInput2->Text = L"Input file 2:";
			// 
			// _btnInput2
			// 
			this->_btnInput2->Location = System::Drawing::Point(259, 57);
			this->_btnInput2->Name = L"_btnInput2";
			this->_btnInput2->Size = System::Drawing::Size(28, 23);
			this->_btnInput2->TabIndex = 53;
			this->_btnInput2->Text = L"..";
			this->_btnInput2->UseVisualStyleBackColor = true;
			this->_btnInput2->Click += gcnew System::EventHandler(this, &MainWindow::_btnInput2_Click);
			// 
			// _btnInput1
			// 
			this->_btnInput1->Location = System::Drawing::Point(259, 31);
			this->_btnInput1->Name = L"_btnInput1";
			this->_btnInput1->Size = System::Drawing::Size(28, 23);
			this->_btnInput1->TabIndex = 52;
			this->_btnInput1->Text = L"..";
			this->_btnInput1->UseVisualStyleBackColor = true;
			this->_btnInput1->Click += gcnew System::EventHandler(this, &MainWindow::_btnInput1_Click);
			// 
			// _txtInput1
			// 
			this->_txtInput1->Location = System::Drawing::Point(79, 33);
			this->_txtInput1->Name = L"_txtInput1";
			this->_txtInput1->Size = System::Drawing::Size(174, 20);
			this->_txtInput1->TabIndex = 50;
			// 
			// _lblInput1
			// 
			this->_lblInput1->AutoSize = true;
			this->_lblInput1->Location = System::Drawing::Point(14, 36);
			this->_lblInput1->Name = L"_lblInput1";
			this->_lblInput1->Size = System::Drawing::Size(59, 13);
			this->_lblInput1->TabIndex = 49;
			this->_lblInput1->Text = L"Input file 1:";
			// 
			// _txtInput2
			// 
			this->_txtInput2->Location = System::Drawing::Point(79, 59);
			this->_txtInput2->Name = L"_txtInput2";
			this->_txtInput2->Size = System::Drawing::Size(174, 20);
			this->_txtInput2->TabIndex = 48;
			// 
			// _lstConsole
			// 
			this->_lstConsole->BackColor = System::Drawing::Color::Black;
			this->_lstConsole->Font = (gcnew System::Drawing::Font(L"Lucida Console", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->_lstConsole->ForeColor = System::Drawing::Color::White;
			this->_lstConsole->FormattingEnabled = true;
			this->_lstConsole->ItemHeight = 12;
			this->_lstConsole->Location = System::Drawing::Point(32, 64);
			this->_lstConsole->Name = L"_lstConsole";
			this->_lstConsole->Size = System::Drawing::Size(428, 28);
			this->_lstConsole->TabIndex = 1;
			// 
			// _mnuContext
			// 
			this->_mnuContext->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {this->_mnuCCompare, this->_mnuCPlot});
			this->_mnuContext->Name = L"_mnuContext";
			this->_mnuContext->Size = System::Drawing::Size(124, 48);
			// 
			// _mnuCCompare
			// 
			this->_mnuCCompare->Name = L"_mnuCCompare";
			this->_mnuCCompare->Size = System::Drawing::Size(123, 22);
			this->_mnuCCompare->Text = L"&Compare";
			this->_mnuCCompare->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCCompare_Click);
			// 
			// _mnuCPlot
			// 
			this->_mnuCPlot->Name = L"_mnuCPlot";
			this->_mnuCPlot->Size = System::Drawing::Size(123, 22);
			this->_mnuCPlot->Text = L"&Plot";
			this->_mnuCPlot->Click += gcnew System::EventHandler(this, &MainWindow::_mnuCPlot_Click);
			// 
			// MainWindow
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(840, 771);
			this->Controls->Add(this->_split);
			this->Controls->Add(this->_tsToolbar);
			this->Controls->Add(this->_stsStatusbar);
			this->Controls->Add(this->_mnuMainMenu);
			this->Icon = (cli::safe_cast<System::Drawing::Icon^  >(resources->GetObject(L"$this.Icon")));
			this->MainMenuStrip = this->_mnuMainMenu;
			this->Name = L"MainWindow";
			this->Opacity = 0;
			this->Text = L"Molecular::Dynamics";
			this->WindowState = System::Windows::Forms::FormWindowState::Maximized;
			this->FormClosing += gcnew System::Windows::Forms::FormClosingEventHandler(this, &MainWindow::MainWindow_FormClosing);
			this->_mnuMainMenu->ResumeLayout(false);
			this->_mnuMainMenu->PerformLayout();
			this->_stsStatusbar->ResumeLayout(false);
			this->_stsStatusbar->PerformLayout();
			this->_tsToolbar->ResumeLayout(false);
			this->_tsToolbar->PerformLayout();
			this->_split->Panel1->ResumeLayout(false);
			this->_split->Panel2->ResumeLayout(false);
			this->_split->ResumeLayout(false);
			this->_pnlSplit->Panel1->ResumeLayout(false);
			this->_pnlSplit->Panel2->ResumeLayout(false);
			this->_pnlSplit->ResumeLayout(false);
			this->_pnlConfigPanel->ResumeLayout(false);
			this->_pnlConfigPanel->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudAnimTS))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudOutputTimesteps))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudMBParticles))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudTimesteps))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_nudParticles))->EndInit();
			this->_pnlCompare->ResumeLayout(false);
			this->_pnlCompareWindow->ResumeLayout(false);
			this->_pnlCompareWindow->PerformLayout();
			this->_mnuContext->ResumeLayout(false);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
};
}

#endif