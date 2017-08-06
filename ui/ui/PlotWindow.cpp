///
/// Plot window controller implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "stdafx.h"
#include "PlotWindow.h"
#include "..\..\common\constants.h"

namespace ui
{
	vector<float*> posVec; //mixed types
	vector<float*> velVec;
	int numberOfParticles;
	int count;
	int ts;

	PlotWindow::PlotWindow(configuration * config)
	{
		_titleBar = gcnew ui::TitleBar();
		_titleBar->Visible = true;
		_titleBar->Location.X = 0;
		_titleBar->Location.Y = 0;
		_titleBar->Dock = System::Windows::Forms::DockStyle::Top;
		_titleBar->CloseClicked += gcnew EventHandler(this, &PlotWindow::_titleBarEvent);
		_titleBar->SetTitle("Plot");
		InitializeComponent();
		
		this->_config = config;
		_play = false;
		
		_nova = NULL;
		_camera = new NovaFX::NovaCamera::Camera(1.15f, 1.15f, 1.15f);

		_particles = gcnew Dictionary<String^, List<int>^>();

		Visible = false;

		PlotWindow_Resize(this, nullptr);

		_tmrDC->Enabled = true;

		_lastMX = 0;
		_lastMY = 0;
		Controls->Add(_titleBar);
	}

	System::Void PlotWindow::_titleBarEvent(System::Object^  , System::EventArgs^  )
	{
		Visible = false;
		if (Parent)
			Parent->Visible = false;
	}

	PlotWindow::~PlotWindow()
	{
		Destroy();
		delete _camera;

		if (components)
		{
			delete components;
		}
	}

	System::Void PlotWindow::PlotWindow_Resize(System::Object^  , System::EventArgs^ ) 
	{
		if (_nova != NULL)
			_nova->Resize(_pnlDC->Width, _pnlDC->Height);
	}

	void PlotWindow::Init()
	{
		if (_nova != NULL)
		{
			delete _nova;
			_nova = NULL;
		}

		HWND handle((HWND) _pnlDC->Handle.ToPointer());

		_nova = new NovaFX::NovaFX(handle, _pnlDC->Width, _pnlDC->Height, 32);
		_nova->SetCamera(_camera);
		_nova->Initialize(0,0,0,0);		
	}

	void PlotWindow::Destroy()
	{
		if (_nova==NULL)
			return;
		
		delete _nova;
		_nova = NULL;
	}

	System::Void PlotWindow::_tmrDC_Tick(System::Object^  , System::EventArgs^  )
	{
		if (_nova == NULL)
			Init();

		if (!Visible)
		{
			_play = false;
			return;
		}

		if (_play)
		{
			if ((_animBar->Value == _animBar->Maximum && !_chkReverse->Checked)
				|| 
				(_animBar->Value == _animBar->Minimum && _chkReverse->Checked))
			{
				_play = false;
				_btnPlayStop->Text = "Play";
			}
			else
			{
				_animBar->Value += (_chkReverse->Checked ? -1 : 1);
				LoadIndex(_animBar->Value, numberOfParticles);
			}
		}

		if (_nova != NULL)
		{
			if (!_nova->Render())
				return;

			if (GetAsyncKeyState('W'))
				_camera->MoveForward(-1);

			if (GetAsyncKeyState('S'))
				_camera->MoveForward(1);

			if (GetAsyncKeyState('A'))
				_camera->StrafeRight(-1);

			if (GetAsyncKeyState('D'))
				_camera->StrafeRight(1);
		}
	}

	System::Void PlotWindow::_pnlDC_MouseDown(System::Object^  , System::Windows::Forms::MouseEventArgs^  e)
	{
		if (e->Button == System::Windows::Forms::MouseButtons::Left)
			_leftButtonState = true;
	}

	System::Void PlotWindow::_pnlDC_MouseUp(System::Object^  , System::Windows::Forms::MouseEventArgs^  e)
	{
		if (e->Button == System::Windows::Forms::MouseButtons::Left)
			_leftButtonState = false;
	}

	System::Void PlotWindow::_pnlDC_MouseMove(System::Object^  , System::Windows::Forms::MouseEventArgs^  e) 
	{
		int xd = e->X - _lastMX;
		int yd = e->Y - _lastMY;

		if (_leftButtonState)
		{
			_camera->RotateY((float)xd / -10.0f);
			_camera->RotateX((float)yd / -10.0f);
		}

		_lastMX = e->X;
		_lastMY = e->Y;
	}

	//open dialog method - return filename on success or null on failure
	System::String ^ PlotWindow::ChooseFile(System::String ^ title, System::String ^ filter)
	{
		OpenFileDialog^ ofd = gcnew OpenFileDialog();
		ofd->InitialDirectory = Application::StartupPath;
		ofd->Filter = filter;
		ofd->FilterIndex = 1;
		ofd->RestoreDirectory = true;
		ofd->Title = title;

		if (ofd->ShowDialog() == System::Windows::Forms::DialogResult::OK)
			return ofd->FileName;

		return nullptr;
	}

	System::Void PlotWindow::_btnBrowse_Click(System::Object^  , System::EventArgs^  ) 
	{
		System::String ^ filename = ChooseFile("Choose an input file", "Animation files |*.anm");

		if (!String::IsNullOrEmpty(filename))
			LoadAnimationFile(MakeCSTR(filename));
	}

	void PlotWindow::flushVectors()
	{
		for(vector<float*>::iterator it=posVec.begin() ; it < posVec.end(); it++)
			delete [] *it;

		for(vector<float*>::iterator it=velVec.begin() ; it < velVec.end(); it++)
			delete [] *it;

		posVec.clear();
		velVec.clear();
	}

	void PlotWindow::LoadAnimationFile(string filename)
	{
		if (_nova==NULL)
			return;

		_tmrDC->Enabled = false;
		_chkVec->Enabled = false;
		_animBar->Enabled = false;
		_btnPrev->Enabled = false;
		_btnNext->Enabled = false;
		_btnPlayStop->Enabled = false;
		_chkReverse->Enabled = false;

		_nova->RemoveAllSpheres();
		_nova->RemoveAllCylinders();
		flushVectors();

		ifstream anim(filename.c_str(), ifstream::in | ifstream::binary);
		if (anim.fail())
		{
			MSGERR("Unable to open input file!");
			return;
		}

		anim.read((char*) &numberOfParticles, sizeof(int));
		anim.read((char*) &ts, sizeof(int));
		anim.read((char*) &count, sizeof(int));

		for(int i=0; i<count; i++)
		{
			float * posArray;

			try
			{
				posArray = new float[numberOfParticles*4];
			}
			catch(bad_alloc)
			{
				MSGERR("Unable to allocate memory!");

				if (posArray)
					delete [] posArray;

				flushVectors();

				_tmrDC->Enabled = true;

				return;
			}

			for(int n=0; n<numberOfParticles*4; n++)
				anim.read((char*) (posArray+n), sizeof(float));

			posVec.push_back(posArray);
		}

		for(int i=0; i<count; i++)
		{
			float * velocityArray;

			try
			{
				velocityArray = new float[numberOfParticles*3];
			}
			catch(bad_alloc)
			{
				MSGERR("Unable to allocate memory!");

				if (velocityArray)
					delete [] velocityArray;

				flushVectors();

				_tmrDC->Enabled = true;

				return;
			}

			for(int n=0; n<numberOfParticles*3; n++)
				anim.read((char*) (velocityArray+n), sizeof(float));

			velVec.push_back(velocityArray);
		}

		LoadIndex(0, numberOfParticles);

		_animBar->Value = 0;
		_animBar->Maximum = count-1;
		_animBar->Enabled = true;
		_btnPrev->Enabled = true;
		_btnNext->Enabled = true;
		_btnPlayStop->Enabled = true;
		_chkReverse->Enabled = true;
		_chkVec->Enabled = true;
		_tmrDC->Enabled = true;
	}

	void PlotWindow::LoadIndex(int index, int nop)
	{
		if (posVec.empty() || velVec.empty() || index > (int)posVec.size() || index > (int)velVec.size())
			return;

		_nova->RemoveAllSpheres();
		_nova->RemoveAllCylinders();

		float * posArray = posVec[index];
		float * velocityArray = velVec[index];

		for(int ip=0; ip<nop; ip++)
		{
			NovaFX::GLRGBA color(0.0f,0.0f,0.0f,1.0f);
			NovaFX::GLRGBA velColor(1.0f,1.0f,1.0f,1.0f);
			float size = 1;

			switch((int)((NovaFX::Nfloat4*)posArray)[ip].w)
			{
			case N:
				color.b = 1.0f;
				size = radN;
				break;
			case O:
				color.r = 1.0f;
				size = radO;
				break;
			case NG:
				color.g = 1.0f;
				size = radNG;
				break;
			default:
				size =0;
				break;
			}

			NovaFX::NSphere * sphere = new NovaFX::NSphere();
			sphere->pos->x = ((NovaFX::Nfloat4*)posArray)[ip].x;
			sphere->pos->y = ((NovaFX::Nfloat4*)posArray)[ip].y;
			sphere->pos->z = ((NovaFX::Nfloat4*)posArray)[ip].z;
			sphere->pos->w = size;
			sphere->slices = sphere->stacks = 16;
			sphere->color = color;
			_nova->AddSphere(sphere);

		}
	}

	void PlotWindow::LoadSimpleFile(string filename)
	{
		if (_nova==NULL)
			return;

		_tmrDC->Enabled = false;
		_chkVec->Enabled = false;
		_animBar->Enabled = false;
		_btnPrev->Enabled = false;
		_btnNext->Enabled = false;
		_btnPlayStop->Enabled = false;
		_chkReverse->Enabled = false;

		_nova->RemoveAllSpheres();
		_nova->RemoveAllCylinders();
		flushVectors();

		float * posArray = NULL;
		float * velocityArray = NULL;

		try
		{
			posArray = new float[(_config->LennardJonesParticles + _config->ManyBodyParticles) * 4];
			velocityArray = new float[(_config->LennardJonesParticles + _config->ManyBodyParticles) * 3];
		}
		catch(bad_alloc)
		{
			MSGERR("Unable to allocate memory!");

			if (posArray)
				delete [] posArray;

			if (velocityArray)
				delete [] velocityArray;

			_tmrDC->Enabled = true;

			return;
		}

		if (!ReadInput(filename.c_str(),(_config->LennardJonesParticles + _config->ManyBodyParticles), posArray, velocityArray))
		{
			MSGERR("Unable to open input file!");

			delete [] posArray;
			delete [] velocityArray;

			return;
		}

		for(int ip=0; ip<(_config->LennardJonesParticles + _config->ManyBodyParticles); ip++)
		{
			NovaFX::GLRGBA color(0.0f,0.0f,0.0f,1.0f);
			float size = 1;

			switch((int)((NovaFX::Nfloat4*)posArray)[ip].w)
			{
			case N:
				color.b = 1.0f;
				size = radN;
				break;
			case O:
				color.r = 1.0f;
				size = radO;
				break;
			case NG:
				color.g = 1.0f;
				size = radNG;
				break;
			}

			NovaFX::NSphere * sphere = new NovaFX::NSphere();
			sphere->pos->x = ((NovaFX::Nfloat4*)posArray)[ip].x;
			sphere->pos->y = ((NovaFX::Nfloat4*)posArray)[ip].y;
			sphere->pos->z = ((NovaFX::Nfloat4*)posArray)[ip].z;
			sphere->pos->w = size;
			sphere->slices = sphere->stacks = 16;
			sphere->color = color;
			_nova->AddSphere(sphere);

		}
		
		delete [] posArray;
		delete [] velocityArray;

		_tmrDC->Enabled = true;
		_chkVec->Enabled = true;
	}

	System::Void PlotWindow::_animBar_Scroll(System::Object^  , System::EventArgs^  )
	{
		LoadIndex(_animBar->Value, numberOfParticles);
	}

	System::Void PlotWindow::_btnPlayStop_Click(System::Object^  , System::EventArgs^  )
	{
		_play = _btnPlayStop->Text->Equals("Play");
		_btnPlayStop->Text = _btnPlayStop->Text->Equals("Play") ? "Stop" : "Play";
	}

	System::Void PlotWindow::_chkVec_CheckedChanged(System::Object^  , System::EventArgs^ e)
	{
		if (_chkVec->Checked)
			_nova->ShowCylinders();
		else
			_nova->HideCylinders();
	}
}