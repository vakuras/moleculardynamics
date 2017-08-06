///
/// Plot window controller definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#pragma once

#ifndef PLOTH
#define PLOTH

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;
using namespace System::Collections::Generic;

#include "..\NovaFX\NovaFX.h"
#include "TitleBar.h"

namespace ui 
{
	public ref class PlotWindow : public System::Windows::Forms::UserControl
	{
	private:
		NovaFX::NovaFX * _nova;
		NovaFX::NovaCamera::Camera * _camera;
		configuration * _config;
		int _lastMX;
		int _lastMY;
		System::Windows::Forms::PictureBox^  _pnlDC;
		System::Windows::Forms::Timer^  _tmrDC;
		TitleBar ^ _titleBar;
		Dictionary<String^, List<int>^>^ _particles;
		TableLayoutPanel^  _tablePlot;
		TableLayoutPanel^  _tableControls;
		CheckBox^  _chkVec;
		TrackBar^  _animBar;
		Button^  _btnPrev;
		Button^  _btnNext;
		Button^  _btnBrowse;
		System::ComponentModel::IContainer^  components;
		System::Windows::Forms::Button^  _btnPlayStop;
		System::Windows::Forms::CheckBox^  _chkReverse;
		bool _leftButtonState;
		bool _play;

	private:
		System::Void _titleBarEvent(System::Object^  , System::EventArgs^  );
		System::Void PlotWindow_Resize(System::Object^  , System::EventArgs^ );
		System::Void _tmrDC_Tick(System::Object^  , System::EventArgs^  );
		System::Void _pnlDC_MouseMove(System::Object^  , System::Windows::Forms::MouseEventArgs^  );
		System::Void PlotWindow::_btnBrowse_Click(System::Object^  , System::EventArgs^  );
		System::Void _pnlDC_MouseDown(System::Object^  , System::Windows::Forms::MouseEventArgs^  );
		System::Void _pnlDC_MouseUp(System::Object^  , System::Windows::Forms::MouseEventArgs^  );
		System::String ^ PlotWindow::ChooseFile(System::String ^ title, System::String ^ filter);
		System::Void _btnPlayStop_Click(System::Object^  , System::EventArgs^  );
		System::Void _animBar_Scroll(System::Object^  , System::EventArgs^  );
		System::Void _chkVec_CheckedChanged(System::Object^  , System::EventArgs^ );
		void PlotWindow::flushVectors();
		void PlotWindow::LoadIndex(int index, int nop);
		

	public:
		void PlotWindow::LoadSimpleFile(string filename);
		void PlotWindow::LoadAnimationFile(string filename);
		PlotWindow(configuration * );
		void Init();
		void Destroy();

	protected:
		~PlotWindow();

#pragma region Windows Form Designer generated code
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->_pnlDC = (gcnew System::Windows::Forms::PictureBox());
			this->_tmrDC = (gcnew System::Windows::Forms::Timer(this->components));
			this->_tablePlot = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->_tableControls = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->_btnPrev = (gcnew System::Windows::Forms::Button());
			this->_btnNext = (gcnew System::Windows::Forms::Button());
			this->_btnBrowse = (gcnew System::Windows::Forms::Button());
			this->_chkVec = (gcnew System::Windows::Forms::CheckBox());
			this->_animBar = (gcnew System::Windows::Forms::TrackBar());
			this->_btnPlayStop = (gcnew System::Windows::Forms::Button());
			this->_chkReverse = (gcnew System::Windows::Forms::CheckBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_pnlDC))->BeginInit();
			this->_tablePlot->SuspendLayout();
			this->_tableControls->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_animBar))->BeginInit();
			this->SuspendLayout();
			// 
			// _pnlDC
			// 
			this->_pnlDC->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_pnlDC->Location = System::Drawing::Point(3, 3);
			this->_pnlDC->Name = L"_pnlDC";
			this->_pnlDC->Size = System::Drawing::Size(574, 447);
			this->_pnlDC->TabIndex = 17;
			this->_pnlDC->TabStop = false;
			this->_pnlDC->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &PlotWindow::_pnlDC_MouseMove);
			this->_pnlDC->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &PlotWindow::_pnlDC_MouseDown);
			this->_pnlDC->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &PlotWindow::_pnlDC_MouseUp);
			// 
			// _tmrDC
			// 
			this->_tmrDC->Enabled = true;
			this->_tmrDC->Interval = 30;
			this->_tmrDC->Tick += gcnew System::EventHandler(this, &PlotWindow::_tmrDC_Tick);
			// 
			// _tablePlot
			// 
			this->_tablePlot->ColumnCount = 1;
			this->_tablePlot->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 100)));
			this->_tablePlot->Controls->Add(this->_pnlDC, 0, 0);
			this->_tablePlot->Controls->Add(this->_tableControls, 0, 1);
			this->_tablePlot->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_tablePlot->Location = System::Drawing::Point(0, 0);
			this->_tablePlot->Name = L"_tablePlot";
			this->_tablePlot->RowCount = 2;
			this->_tablePlot->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 100)));
			this->_tablePlot->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 40)));
			this->_tablePlot->Size = System::Drawing::Size(580, 493);
			this->_tablePlot->TabIndex = 18;
			// 
			// _tableControls
			// 
			this->_tableControls->ColumnCount = 7;
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				30)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				30)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				50)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				80)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				110)));
			this->_tableControls->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				40)));
			this->_tableControls->Controls->Add(this->_btnPrev, 0, 0);
			this->_tableControls->Controls->Add(this->_btnNext, 1, 0);
			this->_tableControls->Controls->Add(this->_btnBrowse, 6, 0);
			this->_tableControls->Controls->Add(this->_chkVec, 5, 0);
			this->_tableControls->Controls->Add(this->_animBar, 4, 0);
			this->_tableControls->Controls->Add(this->_btnPlayStop, 2, 0);
			this->_tableControls->Controls->Add(this->_chkReverse, 3, 0);
			this->_tableControls->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_tableControls->Location = System::Drawing::Point(3, 456);
			this->_tableControls->Name = L"_tableControls";
			this->_tableControls->RowCount = 1;
			this->_tableControls->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 100)));
			this->_tableControls->Size = System::Drawing::Size(574, 34);
			this->_tableControls->TabIndex = 18;
			// 
			// _btnPrev
			// 
			this->_btnPrev->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_btnPrev->Enabled = false;
			this->_btnPrev->Font = (gcnew System::Drawing::Font(L"Lucida Console", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->_btnPrev->Location = System::Drawing::Point(3, 3);
			this->_btnPrev->Name = L"_btnPrev";
			this->_btnPrev->Size = System::Drawing::Size(24, 28);
			this->_btnPrev->TabIndex = 4;
			this->_btnPrev->Text = L"<";
			this->_btnPrev->UseVisualStyleBackColor = true;
			// 
			// _btnNext
			// 
			this->_btnNext->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_btnNext->Enabled = false;
			this->_btnNext->Font = (gcnew System::Drawing::Font(L"Lucida Console", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->_btnNext->Location = System::Drawing::Point(33, 3);
			this->_btnNext->Name = L"_btnNext";
			this->_btnNext->Size = System::Drawing::Size(24, 28);
			this->_btnNext->TabIndex = 5;
			this->_btnNext->Text = L">";
			this->_btnNext->UseVisualStyleBackColor = true;
			// 
			// _btnBrowse
			// 
			this->_btnBrowse->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_btnBrowse->Font = (gcnew System::Drawing::Font(L"Lucida Console", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->_btnBrowse->Location = System::Drawing::Point(537, 3);
			this->_btnBrowse->Name = L"_btnBrowse";
			this->_btnBrowse->Size = System::Drawing::Size(34, 28);
			this->_btnBrowse->TabIndex = 6;
			this->_btnBrowse->Text = L"..";
			this->_btnBrowse->UseVisualStyleBackColor = true;
			this->_btnBrowse->Click += gcnew System::EventHandler(this, &PlotWindow::_btnBrowse_Click);
			// 
			// _chkVec
			// 
			this->_chkVec->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Left | System::Windows::Forms::AnchorStyles::Right));
			this->_chkVec->AutoSize = true;
			this->_chkVec->Checked = true;
			this->_chkVec->CheckState = System::Windows::Forms::CheckState::Checked;
			this->_chkVec->Enabled = false;
			this->_chkVec->Location = System::Drawing::Point(427, 8);
			this->_chkVec->Name = L"_chkVec";
			this->_chkVec->Size = System::Drawing::Size(104, 17);
			this->_chkVec->TabIndex = 1;
			this->_chkVec->Text = L"Show Velocities";
			this->_chkVec->UseVisualStyleBackColor = true;
			this->_chkVec->CheckedChanged += gcnew System::EventHandler(this, &PlotWindow::_chkVec_CheckedChanged);
			// 
			// _animBar
			// 
			this->_animBar->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_animBar->Enabled = false;
			this->_animBar->Location = System::Drawing::Point(193, 3);
			this->_animBar->Name = L"_animBar";
			this->_animBar->Size = System::Drawing::Size(228, 28);
			this->_animBar->TabIndex = 3;
			this->_animBar->Scroll += gcnew System::EventHandler(this, &PlotWindow::_animBar_Scroll);
			// 
			// _btnPlayStop
			// 
			this->_btnPlayStop->Dock = System::Windows::Forms::DockStyle::Fill;
			this->_btnPlayStop->Enabled = false;
			this->_btnPlayStop->Location = System::Drawing::Point(63, 3);
			this->_btnPlayStop->Name = L"_btnPlayStop";
			this->_btnPlayStop->Size = System::Drawing::Size(44, 28);
			this->_btnPlayStop->TabIndex = 7;
			this->_btnPlayStop->Text = L"Play";
			this->_btnPlayStop->UseVisualStyleBackColor = true;
			this->_btnPlayStop->Click += gcnew System::EventHandler(this, &PlotWindow::_btnPlayStop_Click);
			// 
			// _chkReverse
			// 
			this->_chkReverse->Anchor = System::Windows::Forms::AnchorStyles::Left;
			this->_chkReverse->AutoSize = true;
			this->_chkReverse->Enabled = false;
			this->_chkReverse->Location = System::Drawing::Point(113, 8);
			this->_chkReverse->Name = L"_chkReverse";
			this->_chkReverse->Size = System::Drawing::Size(66, 17);
			this->_chkReverse->TabIndex = 8;
			this->_chkReverse->Text = L"Reverse";
			this->_chkReverse->UseVisualStyleBackColor = true;
			// 
			// PlotWindow
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->Controls->Add(this->_tablePlot);
			this->Name = L"PlotWindow";
			this->Size = System::Drawing::Size(580, 493);
			this->Resize += gcnew System::EventHandler(this, &PlotWindow::PlotWindow_Resize);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_pnlDC))->EndInit();
			this->_tablePlot->ResumeLayout(false);
			this->_tableControls->ResumeLayout(false);
			this->_tableControls->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->_animBar))->EndInit();
			this->ResumeLayout(false);

		}
#pragma endregion
};
}


#endif