/********************************************************************************
** Form generated from reading ui file 'MainWindow.ui'
**
** Created: Mon Aug 10 19:21:22 2009
**      by: Qt User Interface Compiler version 4.5.1
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QScrollArea>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_frmMain
{
public:
    QScrollArea *glWidgetArea;
    QWidget *scrollAreaWidgetContents;

    void setupUi(QWidget *frmMain)
    {
        if (frmMain->objectName().isEmpty())
            frmMain->setObjectName(QString::fromUtf8("frmMain"));
        frmMain->resize(527, 525);
        glWidgetArea = new QScrollArea(frmMain);
        glWidgetArea->setObjectName(QString::fromUtf8("glWidgetArea"));
        glWidgetArea->setGeometry(QRect(9, 8, 511, 511));
        QSizePolicy sizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(glWidgetArea->sizePolicy().hasHeightForWidth());
        glWidgetArea->setSizePolicy(sizePolicy);
        glWidgetArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        glWidgetArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        glWidgetArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 509, 509));
        glWidgetArea->setWidget(scrollAreaWidgetContents);

        retranslateUi(frmMain);

        QMetaObject::connectSlotsByName(frmMain);
    } // setupUi

    void retranslateUi(QWidget *frmMain)
    {
        frmMain->setWindowTitle(QApplication::translate("frmMain", "Computer Graphics LAB FrameWork", 0, QApplication::UnicodeUTF8));
        Q_UNUSED(frmMain);
    } // retranslateUi

};

namespace Ui {
    class frmMain: public Ui_frmMain {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
