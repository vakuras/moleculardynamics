///
/// Results window controller implementation.
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "StdAfx.h"
#include "ResultsWindow.h"
#include <vector>

#define RESULTSIZE 10
#define STEP 300.0f

namespace ui
{
	struct double2
	{
		double x;
		double y;
	};

	vector<double2> values[RESULTSIZE]; 

	static const string titles[] = {"Temperature (T)", 
												  "Kinetic Energy (EK)", 
												  "Potentional Energy (EU)",
												  "Total Energy (E)",
												  "Center Of Mass (R.x)",
												  "Center Of Mass (R.y)",
												  "Center Of Mass (R.z)",
												  "Momentum (p.x)",
												  "Momentum (p.y)",
												  "Momentum (p.z)",};

	DrawingPanel::DrawingPanel(void)
	{
		SetStyle(ControlStyles::UserPaint |
			ControlStyles::AllPaintingInWmPaint |
			ControlStyles::DoubleBuffer, true);
	}

	ResultsWindow::ResultsWindow()
	{
		InitializeComponent();

		Visible = false;
		_tmrDraw->Enabled = true;

		_lastMX = _lastMY = 0;
		_zoom = 1;
		_movementx = 0;
		_movementy = 0;
		Controls->Add(_titleBar);
	}

	void ResultsWindow::InitializeComponent()
	{
		components = (gcnew System::ComponentModel::Container());
		System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(ResultsWindow::typeid));
		_pnlView = (gcnew DrawingPanel());
		_tmrDraw = (gcnew Timer(components));
		_lblSize = (gcnew Label());
		_titleBar = gcnew ui::TitleBar();
		_titleBar->Visible = true;
		_titleBar->Location.X = 0;
		_titleBar->Location.Y = 0;
		_titleBar->Dock = DockStyle::Top;
		_titleBar->CloseClicked += gcnew EventHandler(this, &ResultsWindow::_titleBarEvent);
		_titleBar->SetTitle("Results");
		SuspendLayout();
		_pnlView->Location = Drawing::Point(31, 30);
		_pnlView->Name = L"_pnlView";
		_pnlView->Size = Drawing::Size(640, 480);
		_pnlView->TabIndex = 0;
		_pnlView->Paint += gcnew PaintEventHandler(this, &ResultsWindow::_pnlView_Paint);
		_pnlView->MouseMove += gcnew MouseEventHandler(this, &ResultsWindow::_pnlView_MouseMove);
		_tmrDraw->Interval = 33;
		_tmrDraw->Tick += gcnew EventHandler(this, &ResultsWindow::_tmrDraw_Tick);
		_lblSize->Location = Drawing::Point(31, 513);
		_lblSize->Name = L"_lblSize";
		_lblSize->Size = Drawing::Size(640, 18);
		_lblSize->TabIndex = 4;
		_lblSize->TextAlign = Drawing::ContentAlignment::MiddleCenter;
		AutoScaleDimensions = Drawing::SizeF(6, 13);
		AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
		BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
		Controls->Add(_pnlView);
		Controls->Add(_lblSize);
		Name = L"ResultsWindow";
		Size = Drawing::Size(702, 538);
		Resize += gcnew EventHandler(this, &ResultsWindow::ResultsWindow_Resize);
		ResumeLayout(false);
	}

	ResultsWindow::~ResultsWindow()
	{
		if (components)
		{
			delete components;
		}
	}

	Void ResultsWindow::ResultsWindow_Resize(Object^  sender, EventArgs^  e)
	{
		_pnlView->Left = 12;
		_pnlView->Top = _titleBar->Height + 6;
		_pnlView->Width = Width - 24;
		_pnlView->Height = Height - 18 - _lblSize->Height - _titleBar->Height;

		_lblSize->Top = _pnlView->Top + _pnlView->Height + 6;
		_lblSize->Left = Width / 2 - _lblSize->Width / 2;
	}

	Void ResultsWindow::_pnlView_Paint(Object^  , PaintEventArgs^ e)
	{
		try
		{
			Graphics ^ g = e->Graphics;

			g->Clear(Color::Black);

			Pen ^ pen = gcnew Pen(Brushes::White,2.0f);
			Drawing::Font ^ font = gcnew Drawing::Font("Arial", 10);

			for(int i=0; i<RESULTSIZE+1; i++)
				g->DrawLine(pen, 0.0f, i*STEP + _movementy, (float)_pnlView->Width, i*STEP + _movementy);

			for(int i=0; i<RESULTSIZE; i++)
				g->DrawString(gcnew String(titles[i].c_str()), font, Brushes::White, 0, i*STEP + _movementy);

			if (values[0].size()>0)
			{
				//0 = t
				//1 = e
				//2 = ek
				//3 = eu
				//4 = rx
				//5 = ry
				//6 = rz
				//7 = px
				//8 = py
				//9 = pz
				
				double sv[RESULTSIZE];
				double hv[RESULTSIZE];

				array<Pen^>^ pens = gcnew array<Pen^>(RESULTSIZE);

				pens[0] = Pens::Chartreuse;
				pens[1] = Pens::Chocolate;
				pens[2] = Pens::Coral;
				pens[3] = Pens::CornflowerBlue;
				pens[4] = Pens::Cyan;
				pens[5] = Pens::WhiteSmoke;
				pens[6] = Pens::Gold;
				pens[7] = Pens::DodgerBlue;
				pens[8] = Pens::White;
				pens[9] = Pens::Gray;

				for(int i=0; i<RESULTSIZE; i++)
				{
					sv[i] = -1e308;
					hv[i] = 1e308;
				}


				for(unsigned int i=0; i<values[0].size(); i++)
				{
					for(int v=0; v<RESULTSIZE; v++)
					{
						if (values[v][i].y>sv[v])
							sv[v] = values[v][i].y;

						if (values[v][i].y<hv[v])
							hv[v] = values[v][i].y;
					}
				}

				const int val = (int)STEP; 

				for(int v=0; v<RESULTSIZE; v++)
				{
					List<PointF>^ points = gcnew List<PointF>();

					int vxcmCount = -1;

					for(unsigned int i=0; i<values[0].size(); i++)
					{
						PointF point;
						point.X = ((float) values[v][i].x)/_zoom + _movementx;
						double dev = (hv[v]-sv[v] == 0) ? 0.00001 : hv[v]-sv[v];
						point.Y = (float) ((((val * v) + 1) + ((values[v][i].y-sv[v])*(val - 3))/(dev)) + _movementy);
						points->Add(point);
					}

					g->DrawLines(pens[v], points->ToArray()); 

					SizeF ^ size = g->MeasureString(sv[v] + "", font);
					g->DrawString(sv[v] + "", font, Brushes::White, _pnlView->Width - size->Width, (val * v) + _movementy);

					size = g->MeasureString(hv[v] + "", font);
					g->DrawString(hv[v] + "", font, Brushes::White, _pnlView->Width - size->Width, val - size->Height + (val * v) + _movementy);
				}	
			}
		}
		catch(System::Exception^)
		{}
	}

	Void ResultsWindow::_tmrDraw_Tick(Object^  , EventArgs^  )
	{
		if (!Visible)
			return;

		_pnlView->Invalidate(Drawing::Rectangle(0,0,_pnlView->Width, _pnlView->Height));

		_lblSize->Text = "Displacement: " + -_movementx + " Size: 1:" + _zoom;
		_lblSize->Left = Width / 2 - _lblSize->Width / 2;
	}

	bool ResultsWindow::GetResults(String ^ filename)
	{
		_movementx = 0;
		_movementy = 0;
		_zoom=1;
		_lblSize->Text = "Displacement: " + -_movementx + " Size: 1:" + _zoom;
		_lblSize->Left = Width / 2 - _lblSize->Width / 2;

		for(int i=0; i<RESULTSIZE; i++)
			values[i].clear();
		
		ifstream result(MakeCSTR(filename).c_str(), ifstream::in | ifstream::binary);
		int itmp;
		float ftmp;
		unsigned char usebyte;
		double2 tmpvec;

		CPPSAFE_CALL(result.fail(), "Error opening result file!");

		result.read((char*) &itmp, sizeof(int));
		_titleBar->SetTitle("Results [Particles: " + itmp);

		result.read((char*) &itmp, sizeof(int));
		_titleBar->AppendTitle("; Timesteps: " + itmp);

		result.read((char*) &ftmp, sizeof(float));
		_titleBar->AppendTitle(L"; δt: " + ftmp);

		result.read((char*)&usebyte, sizeof(unsigned char));

		if (usebyte == 0xFF)
		{
			result.read((char*) &ftmp, sizeof(float));
			_titleBar->AppendTitle("; LJrs: " + ftmp);

			result.read((char*) &ftmp, sizeof(float));
			_titleBar->AppendTitle("; LJrcut: " + ftmp);
		}

		result.read((char*)&usebyte, sizeof(unsigned char));

		if (usebyte == 0xFF)
		{
			result.read((char*) &ftmp, sizeof(float));
			_titleBar->AppendTitle("; MBrs: " + ftmp);

			result.read((char*) &ftmp, sizeof(float));
			_titleBar->AppendTitle("; MBrcut: " + ftmp);
		}

		result.read((char*)&ftmp, sizeof(float));

		_titleBar->AppendTitle("; vxcm: " + ftmp + ";]");

		//check for errors after read
		CPPSAFE_CALL(result.fail(), "Error reading from result file.");

		for(;;)
		{
			result.read((char*) &tmpvec.x, sizeof(double));

			if (result.eof())
				break;

			//0 = t
			//1 = e
			//2 = ek
			//3 = eu
			//4 = rx
			//5 = ry
			//6 = rz
			//7 = px
			//8 = py
			//9 = pz

			for(int i=0; i<RESULTSIZE; i++)
			{
				result.read((char*) &tmpvec.y, sizeof(double));
				values[i].push_back(tmpvec);
			}

			//check for errors after read
			CPPSAFE_CALL(result.fail(), "Error reading from result file.");
		}

		result.close();
		
		return true;
	}

	Void ResultsWindow::_pnlView_MouseMove(Object^  , MouseEventArgs^  e)
	{
		int xd = e->X - _lastMX;
		int yd = e->Y - _lastMY;

		if (e->Button == System::Windows::Forms::MouseButtons::Left)
		{
			_movementx += (float) xd;
			_movementy += (float) yd;
		}

		if (e->Button == System::Windows::Forms::MouseButtons::Right)
		{
			_zoom += ((float) yd) / 50.0f;

			if (_zoom <= 0)
				_zoom = 0.01f;
		}

		_lastMX = e->X;
		_lastMY = e->Y;
	}

	Void ResultsWindow::_titleBarEvent(Object^  , EventArgs^ )
	{
		Visible = false;
		if (Parent)
			Parent->Visible = false;
	}
}
