/****************************************************************************
** HelpWindow meta object code from reading C++ file 'helpwindow.h'
**
** Created: Thu Jan 24 13:59:38 2002
**      by: The Qt MOC ($Id: moc_helpwindow.cpp,v 1.1 2003-06-04 10:19:35 childcvs Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 9
#elif Q_MOC_OUTPUT_REVISION != 9
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "helpwindow.h"
#include <qmetaobject.h>
#include <qapplication.h>



const char *HelpWindow::className() const
{
    return "HelpWindow";
}

QMetaObject *HelpWindow::metaObj = 0;

void HelpWindow::initMetaObject()
{
    if ( metaObj )
	return;
    if ( qstrcmp(QMainWindow::className(), "QMainWindow") != 0 )
	badSuperclassWarning("HelpWindow","QMainWindow");
    (void) staticMetaObject();
}

#ifndef QT_NO_TRANSLATION

QString HelpWindow::tr(const char* s)
{
    return qApp->translate( "HelpWindow", s, 0 );
}

QString HelpWindow::tr(const char* s, const char * c)
{
    return qApp->translate( "HelpWindow", s, c );
}

#endif // QT_NO_TRANSLATION

QMetaObject* HelpWindow::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    (void) QMainWindow::staticMetaObject();
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    typedef void (HelpWindow::*m1_t0)(bool);
    typedef void (QObject::*om1_t0)(bool);
    typedef void (HelpWindow::*m1_t1)(bool);
    typedef void (QObject::*om1_t1)(bool);
    typedef void (HelpWindow::*m1_t2)();
    typedef void (QObject::*om1_t2)();
    typedef void (HelpWindow::*m1_t3)();
    typedef void (QObject::*om1_t3)();
    typedef void (HelpWindow::*m1_t4)();
    typedef void (QObject::*om1_t4)();
    typedef void (HelpWindow::*m1_t5)();
    typedef void (QObject::*om1_t5)();
    typedef void (HelpWindow::*m1_t6)();
    typedef void (QObject::*om1_t6)();
    typedef void (HelpWindow::*m1_t7)();
    typedef void (QObject::*om1_t7)();
    typedef void (HelpWindow::*m1_t8)(const QString&);
    typedef void (QObject::*om1_t8)(const QString&);
    typedef void (HelpWindow::*m1_t9)(int);
    typedef void (QObject::*om1_t9)(int);
    typedef void (HelpWindow::*m1_t10)(int);
    typedef void (QObject::*om1_t10)(int);
    typedef void (HelpWindow::*m1_t11)();
    typedef void (QObject::*om1_t11)();
    m1_t0 v1_0 = &HelpWindow::setBackwardAvailable;
    om1_t0 ov1_0 = (om1_t0)v1_0;
    m1_t1 v1_1 = &HelpWindow::setForwardAvailable;
    om1_t1 ov1_1 = (om1_t1)v1_1;
    m1_t2 v1_2 = &HelpWindow::textChanged;
    om1_t2 ov1_2 = (om1_t2)v1_2;
    m1_t3 v1_3 = &HelpWindow::about;
    om1_t3 ov1_3 = (om1_t3)v1_3;
    m1_t4 v1_4 = &HelpWindow::aboutQt;
    om1_t4 ov1_4 = (om1_t4)v1_4;
    m1_t5 v1_5 = &HelpWindow::openFile;
    om1_t5 ov1_5 = (om1_t5)v1_5;
    m1_t6 v1_6 = &HelpWindow::newWindow;
    om1_t6 ov1_6 = (om1_t6)v1_6;
    m1_t7 v1_7 = &HelpWindow::print;
    om1_t7 ov1_7 = (om1_t7)v1_7;
    m1_t8 v1_8 = &HelpWindow::pathSelected;
    om1_t8 ov1_8 = (om1_t8)v1_8;
    m1_t9 v1_9 = &HelpWindow::histChosen;
    om1_t9 ov1_9 = (om1_t9)v1_9;
    m1_t10 v1_10 = &HelpWindow::bookmChosen;
    om1_t10 ov1_10 = (om1_t10)v1_10;
    m1_t11 v1_11 = &HelpWindow::addBookmark;
    om1_t11 ov1_11 = (om1_t11)v1_11;
    QMetaData *slot_tbl = QMetaObject::new_metadata(12);
    QMetaData::Access *slot_tbl_access = QMetaObject::new_metaaccess(12);
    slot_tbl[0].name = "setBackwardAvailable(bool)";
    slot_tbl[0].ptr = (QMember)ov1_0;
    slot_tbl_access[0] = QMetaData::Private;
    slot_tbl[1].name = "setForwardAvailable(bool)";
    slot_tbl[1].ptr = (QMember)ov1_1;
    slot_tbl_access[1] = QMetaData::Private;
    slot_tbl[2].name = "textChanged()";
    slot_tbl[2].ptr = (QMember)ov1_2;
    slot_tbl_access[2] = QMetaData::Private;
    slot_tbl[3].name = "about()";
    slot_tbl[3].ptr = (QMember)ov1_3;
    slot_tbl_access[3] = QMetaData::Private;
    slot_tbl[4].name = "aboutQt()";
    slot_tbl[4].ptr = (QMember)ov1_4;
    slot_tbl_access[4] = QMetaData::Private;
    slot_tbl[5].name = "openFile()";
    slot_tbl[5].ptr = (QMember)ov1_5;
    slot_tbl_access[5] = QMetaData::Private;
    slot_tbl[6].name = "newWindow()";
    slot_tbl[6].ptr = (QMember)ov1_6;
    slot_tbl_access[6] = QMetaData::Private;
    slot_tbl[7].name = "print()";
    slot_tbl[7].ptr = (QMember)ov1_7;
    slot_tbl_access[7] = QMetaData::Private;
    slot_tbl[8].name = "pathSelected(const QString&)";
    slot_tbl[8].ptr = (QMember)ov1_8;
    slot_tbl_access[8] = QMetaData::Private;
    slot_tbl[9].name = "histChosen(int)";
    slot_tbl[9].ptr = (QMember)ov1_9;
    slot_tbl_access[9] = QMetaData::Private;
    slot_tbl[10].name = "bookmChosen(int)";
    slot_tbl[10].ptr = (QMember)ov1_10;
    slot_tbl_access[10] = QMetaData::Private;
    slot_tbl[11].name = "addBookmark()";
    slot_tbl[11].ptr = (QMember)ov1_11;
    slot_tbl_access[11] = QMetaData::Private;
    metaObj = QMetaObject::new_metaobject(
	"HelpWindow", "QMainWindow",
	slot_tbl, 12,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    metaObj->set_slot_access( slot_tbl_access );
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    return metaObj;
}
