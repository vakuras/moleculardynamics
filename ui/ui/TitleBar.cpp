///
/// Title Bar implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "stdafx.h"
#include "TitleBar.h"

namespace ui 
{
	TitleBar::TitleBar(void)
	{
		resources = (gcnew System::ComponentModel::ComponentResourceManager(TitleBar::typeid));
		InitializeComponent();
	}

	void TitleBar::InitializeComponent()
	{
		_picClose = (gcnew PictureBox());
		_picSave = (gcnew PictureBox());
		_lblTitle = (gcnew Label());
		(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(_picClose))->BeginInit();
		(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(_picSave))->BeginInit();
		SuspendLayout();
		_picClose->BackColor = Drawing::Color::Transparent;
		_picClose->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"_picClose.BackgroundImage")));
		_picClose->BackgroundImageLayout = ImageLayout::None;
		_picClose->Location = Drawing::Point(284, 4);
		_picClose->Name = L"_picClose";
		_picClose->Size = Drawing::Size(31, 17);
		_picClose->TabIndex = 0;
		_picClose->TabStop = false;
		_picClose->MouseLeave += gcnew EventHandler(this, &TitleBar::_picClose_MouseLeave);
		_picClose->MouseDown += gcnew MouseEventHandler(this, &TitleBar::_picClose_MouseDown);
		_picClose->MouseUp += gcnew MouseEventHandler(this, &TitleBar::_picClose_MouseUp);
		_picClose->MouseEnter += gcnew EventHandler(this, &TitleBar::_picClose_MouseEnter);
		_picClose->Click += gcnew EventHandler(this, &TitleBar::_picClose_Click);
		_picSave->BackColor = Drawing::Color::Transparent;
		_picSave->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"save")));
		_picSave->BackgroundImageLayout = ImageLayout::None;
		_picSave->Location = Drawing::Point(284, 4);
		_picSave->Name = L"_picSave";
		_picSave->Size = Drawing::Size(31, 17);
		_picSave->TabIndex = 0;
		_picSave->TabStop = false;
		_picSave->MouseLeave += gcnew EventHandler(this, &TitleBar::_picSave_MouseLeave);
		_picSave->MouseDown += gcnew MouseEventHandler(this, &TitleBar::_picSave_MouseDown);
		_picSave->MouseUp += gcnew MouseEventHandler(this, &TitleBar::_picSave_MouseUp);
		_picSave->MouseEnter += gcnew EventHandler(this, &TitleBar::_picSave_MouseEnter);
		_picSave->Click += gcnew EventHandler(this, &TitleBar::_picSave_Click);
		_picSave->Visible = false;
		_lblTitle->AutoSize = true;
		_lblTitle->BackColor = Drawing::Color::Transparent;
		_lblTitle->Font = (gcnew Drawing::Font(L"Lucida Console", 8.25F, static_cast<Drawing::FontStyle>((Drawing::FontStyle::Bold | Drawing::FontStyle::Italic)), 
			Drawing::GraphicsUnit::Point, static_cast<Byte>(0)));
		_lblTitle->ForeColor = Drawing::Color::Black;
		_lblTitle->Location = Drawing::Point(4, 6);
		_lblTitle->Name = L"_lblTitle";
		_lblTitle->Size = Drawing::Size(37, 11);
		_lblTitle->TabIndex = 1;
		_lblTitle->Text = L"No Title";
		AutoScaleDimensions = Drawing::SizeF(6, 13);
		AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
		BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
		Controls->Add(_lblTitle);
		Controls->Add(_picClose);
		Controls->Add(_picSave);
		MaximumSize = Drawing::Size(200000, 24);
		MinimumSize = Drawing::Size(0, 24);
		Name = L"TitleBar";
		Size = Drawing::Size(318, 24);
		Resize += gcnew EventHandler(this, &TitleBar::TitleBar_Resize);
		(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(_picClose))->EndInit();
		(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(_picSave))->EndInit();
		ResumeLayout(false);
		PerformLayout();
	}

	TitleBar::~TitleBar()
	{
		if (components)
		{
			delete components;
		}
	}

	Void TitleBar::TitleBar_Resize(Object^  sender, EventArgs^  e) 
	{
		_picClose->Left = Width - 37;
		_picSave->Left = _picClose->Left - 6 - _picSave->Width;
	}

	Void TitleBar::_picClose_MouseDown(Object^  sender, MouseEventArgs^  e)
	{
		_picClose->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"close_click")));
	}

	Void TitleBar::_picClose_MouseUp(Object^  sender, MouseEventArgs^  e) 
	{
		_picClose->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"_picClose.BackgroundImage")));
	}

	Void TitleBar::_picClose_MouseEnter(Object^  sender, EventArgs^  e)
	{
		_picClose->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"close_hover")));
	}

	Void TitleBar::_picClose_MouseLeave(Object^  sender, EventArgs^  e)
	{
		_picClose->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"_picClose.BackgroundImage")));
	}

	Void TitleBar::_picClose_Click(Object^  sender, EventArgs^  e)
	{
		CloseClicked(sender, e);
	}

	Void TitleBar::_picSave_MouseDown(Object^  sender, MouseEventArgs^  e)
	{
		_picSave->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"save_click")));
	}

	Void TitleBar::_picSave_MouseUp(Object^  sender, MouseEventArgs^  e) 
	{
		_picSave->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"save")));
	}

	Void TitleBar::_picSave_MouseEnter(Object^  sender, EventArgs^  e)
	{
		_picSave->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"save_hover")));
	}

	Void TitleBar::_picSave_MouseLeave(Object^  sender, EventArgs^  e)
	{
		_picSave->BackgroundImage = (cli::safe_cast<Drawing::Image^  >(resources->GetObject(L"save")));
	}

	Void TitleBar::_picSave_Click(Object^  sender, EventArgs^  e)
	{
		SaveClicked(sender, e);
	}

	Void TitleBar::SetTitle(String ^ title)
	{
		_lblTitle->Text = title;
	}

	Void TitleBar::AppendTitle(String ^ title)
	{
		_lblTitle->Text += title;
	}

	Void TitleBar::SaveSetVisible(bool set)
	{
		_picSave->Visible = set;
	}
}