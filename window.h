
#ifndef WINDOW_H
#define WINDOW_H

#include <QDialog>
#include "ui_MainWindow.h"

class GLWidget;

//Create a window and load the components from generated code ui_MainWindow.h
class Window : public QDialog, public Ui::frmMain
{
	Q_OBJECT;

public:
	//Constructor 
	Window(QWidget *parent = 0);
	
private:
	//Member variable to hold the GLWidget instance
        //We need a reference to this
	GLWidget *m_glWidget;		
};


#endif
