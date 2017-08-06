///
/// Logo window implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "stdafx.h"
#include "LogoWindow.h"
#include "MainWindow.h"

namespace ui
{
	///
	/// Retrieve and create an _instance of the form
	/// 
	LogoWindow ^ LogoWindow::GetInstance(Form ^ owner)
	{
		if (_instance==nullptr)
		{
			_instance = gcnew LogoWindow();
			_instance->Owner = owner;
		}

		return _instance;
	}

	///
	/// Constructor
	///
	LogoWindow::LogoWindow()
	{
		System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(LogoWindow::typeid));
		_tmr = gcnew System::Windows::Forms::Timer();
		_tmr->Interval = 10;
		_tmr->Tick += gcnew System::EventHandler(this, &LogoWindow::_tmr_Tick);
		this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"logo")));
		this->SuspendLayout();
		this->FormBorderStyle = System::Windows::Forms::FormBorderStyle::None;
		this->Opacity = 0.75;
		this->Name = L"LogoWindow";
		this->StartPosition = System::Windows::Forms::FormStartPosition::CenterScreen;
		this->ShowInTaskbar = false;
		this->TopMost = true;
		this->ResumeLayout(false);

		//visuals
		this->BackgroundImageLayout = System::Windows::Forms::ImageLayout::None;
		this->Width = this->BackgroundImage->Width;
		this->Height = this->BackgroundImage->Height;
		this->TransparencyKey = Color::Magenta;

		this->Show();

		percent = 0;
		_tmr->Enabled = true;
	}

	///
	/// Destructor
	///
	LogoWindow::~LogoWindow()
	{
		//dispose of _instance
		_instance = nullptr;

		if (_components)
		{
			delete _components;
		}
	}

	///
	/// Timer tick event
	///
	Void LogoWindow::_tmr_Tick(Object^  , EventArgs^  ) 
	{
		percent+=0.01;
		((MainWindow ^) Owner)->Opacity = percent;

		if (percent>2)
		{
			_tmr->Enabled = false;
			this->Close();
		}
	}
}