///
/// About window definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef ABOUTWIN
#define ABOUTWIN

#pragma once

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;

#include "TitleBar.h"

namespace ui 
{
	public ref class AboutWindow : public Form
	{
	public:
		AboutWindow();

	protected:
		~AboutWindow();
		void InitializeComponent();

	private: 
		System::ComponentModel::Container ^components;
		PictureBox^  _picAuthor;
		Button^  _btnOk;
		Label^  _lblText;
		TitleBar ^ _titleBar;
		bool drag;
		Panel^  _pnlWin;
		Point start;

	private:
		Void _titleBar_MouseDown(Object^  , MouseEventArgs^ e );
		Void _titleBar_MouseMove(Object^  , MouseEventArgs^  e);
		Void _titleBar_MouseUp(Object^  , MouseEventArgs^  );
		Void _titleBarEvent(Object^ , EventArgs^  );
};
}

#endif