
#include <QtGui>
#include <QtOpenGL>

#include "glwidget.h"
#include "window.h"


//------------------------------------------------------------------------------------
// Creates and initializes the main window for application
//------------------------------------------------------------------------------------
Window::Window(QWidget *parent):QDialog(parent)
{
        //We create an instance of GLWidget component we built in glwidget.h
	m_glWidget = new GLWidget;

	//Setup application interface. Creates all the required components and sliders.
	setupUi(this);

	    //We need to attach our m_glWidget to glWidgetArea
        //All our drawings will be on glWidgetArea
        glWidgetArea->setWidget(m_glWidget);
}
