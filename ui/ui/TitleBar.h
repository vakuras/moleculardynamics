///
/// Title Bar definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#ifndef TITLEBAR_H
#define TITLEBAR_H

#pragma once

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;

namespace ui 
{
	public ref class TitleBar : public UserControl
	{
	public:
		event EventHandler^ CloseClicked;
		event EventHandler^ SaveClicked;
		Void SetTitle(String ^ title);
		Void AppendTitle(String ^ title);
		Void SaveSetVisible(bool set);
		TitleBar();

	protected:
		~TitleBar();
		void InitializeComponent();

	private:
		Void TitleBar_Resize(Object^  sender, EventArgs^  e) ;
		Void _picClose_MouseDown(Object^  sender, MouseEventArgs^  e);
		Void _picClose_MouseUp(Object^  sender, MouseEventArgs^  e) ;
		Void _picClose_MouseEnter(Object^  sender, EventArgs^  e);
		Void _picClose_MouseLeave(Object^  sender, EventArgs^  e);
		Void _picClose_Click(Object^  sender, EventArgs^  e);
		Void _picSave_MouseDown(Object^  sender, MouseEventArgs^  e);
		Void _picSave_MouseUp(Object^  sender, MouseEventArgs^  e) ;
		Void _picSave_MouseEnter(Object^  sender, EventArgs^  e);
		Void _picSave_MouseLeave(Object^  sender, EventArgs^  e);
		Void _picSave_Click(Object^  sender, EventArgs^  e);
		PictureBox^  _picClose;
		PictureBox^  _picSave;
		Label^  _lblTitle;
		System::ComponentModel::ComponentResourceManager^  resources;
		System::ComponentModel::Container ^components;
	};
}

#endif