///
/// Logo window definition.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#pragma once

#ifndef LOGOH
#define LOGOH

using namespace System;
using namespace System::Windows::Forms;
using namespace System::Drawing;

namespace ui 
{
	//clr logo window - singleton
	public ref class LogoWindow : public Form
	{
	private:
		//hold the _instance of the form
		static LogoWindow ^ _instance = nullptr;
		System::ComponentModel::Container ^_components;
		System::Windows::Forms::Timer ^ _tmr; //timer for form invalidation
		Void _tmr_Tick(System::Object^  , System::EventArgs^  );
		double percent;

	public:
		static LogoWindow ^ GetInstance(Form ^); //returns instance of the form
		
	protected:
		LogoWindow();
		~LogoWindow();

	};
}

#endif