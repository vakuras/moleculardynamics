///
/// Result window controller definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#pragma once

#ifndef RESULTH
#define RESULTH

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;
using namespace System::Collections::Generic;

#include "TitleBar.h"

namespace ui 
{
	public ref class DrawingPanel : public Panel
	{
	public:
		DrawingPanel();
	};

	public ref class ResultsWindow : public UserControl
	{
	public:
		ResultsWindow();
		bool GetResults(String ^);

	protected:
		~ResultsWindow();
		void InitializeComponent();

	private:
		Void _tmrDraw_Tick(Object^  , EventArgs^  );
		Void _pnlView_MouseMove(Object^  , MouseEventArgs^  );
		Void _pnlView_Paint(Object^  , PaintEventArgs ^  );
		Void ResultsWindow_Resize(Object^  sender, EventArgs^  e);
		Void _titleBarEvent(Object^  , EventArgs^ );

	private: 
		DrawingPanel^  _pnlView;
		Timer^  _tmrDraw;
		Label^  _lblSize;
		System::ComponentModel::IContainer^  components;
		float _movementx;
		float _movementy;
		float _zoom;
		int _lastMX;
		int _lastMY;
		TitleBar ^ _titleBar;
	};
}

#endif