
#include "chintz_base.h"

#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qframe.h>
#include <qgroupbox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qradiobutton.h>
#include <qspinbox.h>
#include <qtabwidget.h>
#include <qtoolbutton.h>
#include <qwidget.h>
#include <qlayout.h>
#include <qvariant.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <qimage.h>
#include <qpixmap.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include <qfileinfo.h>
#include <qfiledialog.h>
#include <string>
#include "helpwindow.h"


static const char* const image0_data[] = { 
"16 13 5 1",
". c None",
"# c #000000",
"b c #353535",
"c c #717171",
"a c #ffffff",
"................",
"...#####........",
"...#aaaa#b......",
"...#aaaa#cb.....",
"...#aaaa#ccb....",
"...#aaaa#####...",
"...#aaaaaaaa#...",
"...#aaaaaaaa#...",
"...#aaaaaaaa#...",
"...#aaaaaaaa#...",
"...#aaaaaaaa#...",
"...##########...",
"................"};

static const char* const image1_data[] = { 
"14 14 4 1",
"b c None",
". c #040404",
"# c #808304",
"a c #bfc2bf",
"..............",
".#.aaaaaaaa.a.",
".#.aaaaaaaa...",
".#.aaaaaaaa.#.",
".#.aaaaaaaa.#.",
".#.aaaaaaaa.#.",
".#.aaaaaaaa.#.",
".##........##.",
".############.",
".##.........#.",
".##......aa.#.",
".##......aa.#.",
".##......aa.#.",
"b............."};

static const char* const image2_data[] = { 
"16 13 5 1",
". c None",
"# c #040404",
"c c #808304",
"a c #f3f704",
"b c #f3f7f3",
".........###....",
"........#...#.#.",
".............##.",
".###........###.",
"#aba#######.....",
"#babababab#.....",
"#ababababa#.....",
"#baba###########",
"#aba#ccccccccc#.",
"#ba#ccccccccc#..",
"#a#ccccccccc#...",
"##ccccccccc#....",
"###########....."};


/* 
 *  Constructs a Chintz_base which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
Chintz_base::Chintz_base( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "Chintz_base" );
    resize( 790, 620 ); 
    QFont f( font() );
    f.setFamily( "adobe-helvetica" );
    setFont( f ); 
    setCaption( tr( "Chintz v0.2"  ) );
    QToolTip::add(  this, tr( "" ) );

    altered_state=false;

    //Make the menu bar at the top of the widget

    menu=new QMenuBar(this);
    menu->setSeparator(QMenuBar::InWindowsStyle);
    QPopupMenu *file=new QPopupMenu(this);
    file->insertItem("&New Control File",this,SLOT(Clear()),CTRL+Key_N);
    file->insertItem("&Open Control File",this,SLOT(OpenChildControlFile()),CTRL+Key_O);
    file->insertItem("&Save Control File",this,SLOT(WriteTheOutput()),CTRL+Key_S);
    file->insertItem("&Save Control File As...",this,SLOT(SaveAs()));
    file->insertItem("&Quit",this,SLOT(Close()),CTRL+Key_Q);
    QPopupMenu *help=new QPopupMenu(this);
    menu->insertSeparator(0);
    help->insertItem("&About",this,SLOT(About()));
    help->insertItem("&Contents",this,SLOT(Help()));
    menu->insertItem("&File",file,4); 
    menu->insertItem("&Help",help,7);
 
    QPixmap image0( ( const char** ) image0_data );
    QPixmap image1( ( const char** ) image1_data );
    QPixmap image2( ( const char** ) image2_data );

    Framebutt = new QFrame( this, "Framebutt" );
    Framebutt->setGeometry( QRect( 10, 30, 750, 40 ) ); 
    Framebutt->setFrameShape( QFrame::Panel );
    Framebutt->setFrameShadow( QFrame::Raised );

    Newby = new QPushButton( Framebutt, "Newby" );
    Newby->setGeometry( QRect( 10, 10, 24, 22 ) ); 
    Newby->setText( tr( ""  ) );
    Newby->setPixmap( image0 );
    QToolTip::add(  Newby, tr( "New Control File" ) );

    Save = new QPushButton( Framebutt, "Save" );
    Save->setGeometry( QRect( 70, 10, 22, 22 ) ); 
    Save->setText( tr( ""  ) );
    Save->setPixmap( image1 );
    QToolTip::add(  Save, tr( "Save Control File" ) );

    OPeny = new QPushButton( Framebutt, "OPeny" );
    OPeny->setGeometry( QRect( 40, 10, 24, 22 ) ); 
    OPeny->setText( tr( ""  ) );
    OPeny->setPixmap( image2 );
    QToolTip::add(  OPeny, tr( "Open Control File" ) );

    connect( Newby, SIGNAL( clicked() ), this, SLOT( Clear() ));
    connect( Save, SIGNAL( clicked() ), this, SLOT( WriteTheOutput() ) );
    connect( OPeny, SIGNAL( clicked() ), this, SLOT( OpenChildControlFile() ));

    CONTROLFILENAME="";
    CONTROLFILENAMELine=new QLineEdit(Framebutt,"CNTRL");
    CONTROLFILENAMELine->setGeometry( QRect( 490, 10, 224, 22 ) ); 
    CONTROLFILENAMELine->setText( tr( ""  ) );
    connect(CONTROLFILENAMELine, SIGNAL( textChanged(const QString&) ), this, SLOT( altered() ) );
    connect(CONTROLFILENAMELine, SIGNAL( returnPressed() ), this, SLOT( control() ) );

    QLabel* Cntrl=new QLabel(Framebutt,"Cntrlabel");
    Cntrl->setGeometry( QRect( 380, 10, 100, 22 ) );
    Cntrl->setText(tr("Child Control File"));
    QFont Cf=Cntrl->font();
    Cf.setBold(true);
    Cntrl->setFont(Cf);

    Br = new QToolButton( Framebutt, "Br" );
    Br->setGeometry( QRect( 714, 10, 20, 23 ) ); 
    Br->setText( tr( "..."  ) );
    QToolTip::add(  Br, tr( "browse" ) );
    connect(Br, SIGNAL( clicked() ), this, SLOT( ChildFileSelect() ) );

    QPalette pal;
    QColorGroup cg;
    cg.setColor( QColorGroup::Foreground, QColor( 200, 110, 250) );
    cg.setColor( QColorGroup::Background, QColor( 192, 192, 192) );
    pal.setActive(cg);
    pal.setInactive(cg);
    Cntrl->setPalette(pal);

    StatusFrame = new QFrame( this, "StatusFrame" );
    StatusFrame->setGeometry( QRect( 10, 590, 751, 21 ) ); 
    StatusFrame->setFrameShape( QFrame::StyledPanel );
    StatusFrame->setFrameShadow( QFrame::Raised );

    Status = new QLabel( StatusFrame, "Status" );
    Status->setGeometry( QRect( 10, 1, 101, 16 ) ); 
    cg.setColor( QColorGroup::Foreground, QColor( 255, 0, 0) );
    cg.setColor( QColorGroup::Background, QColor( 192, 192, 192) );
    pal.setActive(cg);
    pal.setInactive(cg);
    Status->setPalette(pal);

#include "generated_chintz_base_constructor_insert"

    QFileInfo q(".ChintzDefaults");
    if (q.exists()){
      ReadChildControlFile(".ChintzDefaults");
      altered_state=false;
      Status->setText("");
    }
}

/*  
 *  Destroys the object and frees any allocated resources
 */
Chintz_base::~Chintz_base()
{
    // no need to delete child widgets, Qt does it all for us
  ofile.close();
}

/*  
 *  Main event handler. Reimplemented to handle application
 *  font changes
 */


bool Chintz_base::event( QEvent* ev )
{
    bool ret = QDialog::event( ev ); 

    return ret;
}

void Chintz_base::closeEvent( QCloseEvent* ev )
{
  ev->ignore();
  Close();
}


void Chintz_base::Close()
{  if (altered_state){
    QMessageBox *q=new QMessageBox("Really Quit?",
				   "Quit without saving?",
				   QMessageBox::Warning,
				   QMessageBox::Yes,
				   QMessageBox::No,
				   QMessageBox::Cancel,this);
    q->setButtonText( QMessageBox::Yes, "Save" );
    q->setButtonText( QMessageBox::No, "Don't Save" );
    
    switch(q->exec()){
    case QMessageBox::Yes:WriteTheOutput();if (!CONTROLFILENAME.isEmpty())qApp->quit();break;
    case QMessageBox::No:qApp->quit();
    }
  }
  else{
    QMessageBox *q=new QMessageBox("Really Quit?",
				   "Quit?",
				   QMessageBox::Warning,
				   QMessageBox::Yes,
				   QMessageBox::No,
				   QMessageBox::NoButton,this);
    if (q->exec()==QMessageBox::Yes)qApp->quit();
  }
}
void Chintz_base::OpenChildControlFile()
{
  if (altered_state)
    {
    QMessageBox *q=new QMessageBox("Altered state!",
				   "Save the current values?",
				   QMessageBox::Warning,
				   QMessageBox::Yes,
				   QMessageBox::No,
				   QMessageBox::Cancel,this);
    q->setButtonText( QMessageBox::Yes, "Save" );
    q->setButtonText( QMessageBox::No, "Don't Save" );
    
    switch(q->exec()){
    case QMessageBox::Yes:WriteTheOutput();if (CONTROLFILENAME.length() ==0)return;break;
    case QMessageBox::No:break;
    case QMessageBox::Cancel:return;
    }
  }
  QString S=QFileDialog::getOpenFileName(QString::null,"",this);
  if (!S.isEmpty()) {set_CONTROLFILENAME(S);CONTROLFILENAMELine->setText(S);}
  if (!S.isEmpty()){Blank();ReadChildControlFile(CONTROLFILENAME);altered_state=false;Status->setText("");}
}

void Chintz_base::ChildFileSelect(){
  QString S=QFileDialog::getSaveFileName(QString::null,"*",this);
  if (S.length()!=0){set_CONTROLFILENAME(S);CONTROLFILENAMELine->setText(S);}
}

bool Chintz_base::Clear(){
  if (altered_state)
    {
    QMessageBox *q=new QMessageBox("Altered state!",
				   "Save the current values?",
				   QMessageBox::Warning,
				   QMessageBox::Yes,
				   QMessageBox::No,
				   QMessageBox::Cancel,this);
    q->setButtonText( QMessageBox::Yes, "Save" );
    q->setButtonText( QMessageBox::No, "Don't Save" );
    
    switch(q->exec()){
    case QMessageBox::Yes:WriteTheOutput();if (CONTROLFILENAME.length() ==0)return false;break;
    case QMessageBox::No:break;
    case QMessageBox::Cancel:return false;
    }
  }
  CONTROLFILENAME="";CONTROLFILENAMELine->setText("");
  Blank();
  altered_state=false;
  Status->setText("");
  return true;
}

void Chintz_base::SaveAs(){
  QString S=QFileDialog::getSaveFileName(CONTROLFILENAME,"",this);
  if (!S.isEmpty()){
    set_CONTROLFILENAME(S);CONTROLFILENAMELine->setText(S);
  }else{
    return;
  }
  if (CONTROLFILENAME.length() !=0){
    QFileInfo q(CONTROLFILENAME);
    if (q.exists()){
      QString M="File z Exists:Overwite?";
      M.replace(QRegExp("z"),CONTROLFILENAME);
      QMessageBox *q=new QMessageBox("File Exists",
				     M,
				     QMessageBox::Warning,
				     QMessageBox::Yes,
				     QMessageBox::No,
				     QMessageBox::NoButton,this);
      
      switch(q->exec()){
      case QMessageBox::Yes:  ofile.open(CONTROLFILENAME);
	WriteChildControlFile(ofile);
	ofile.close();
	altered_state=false;
	Status->setText("");
      }
    }else{
      ofile.open(CONTROLFILENAME);
      WriteChildControlFile(ofile);
      ofile.close();
      altered_state=false;
      Status->setText("");
    }
  }
}

void Chintz_base::WriteTheOutput()
{
  ofile.close();
  if (CONTROLFILENAME.isEmpty()){
    SaveAs();
  }

  else{
    ofile.open(CONTROLFILENAME);
    WriteChildControlFile(ofile);
    ofile.close();
    altered_state=false;
    Status->setText("");
  }
}

void Chintz_base::About()
{
  QMessageBox *q=new QMessageBox("About Chintz",
			       "Chintz version 0.2\n\nChild Model User Interface\n\nby \n Mike Bithell 2002",
				   QMessageBox::Information,
				   QMessageBox::Ok,
				   QMessageBox::NoButton,
				   QMessageBox::NoButton,this);
 q->exec();
}
void Chintz_base::Help(){
  HelpWindow *help = new HelpWindow("index.html", "docs", 0, "help viewer");
  help->show(); 
}

void Chintz_base::control(){
  QString S=CONTROLFILENAMELine->text();
  QFileInfo q(S);
  QDir d=q.dir();
  if(!d.exists()){
    QString M="Directory z Does not Exist";
    M.replace(QRegExp("z"),d.path());
    QMessageBox *b=new QMessageBox("Error:No Directory",
				   M,
				   QMessageBox::Warning,
				   QMessageBox::Ok,
				   QMessageBox::NoButton,
				   QMessageBox::NoButton,this);
    b->exec();
    CONTROLFILENAMELine->setText(CONTROLFILENAME);  
    return;
  }
  if (q.exists()){
    QString M="File z Exists:What Now?";
    M.replace(QRegExp("z"),CONTROLFILENAME);
    QMessageBox *q=new QMessageBox("File Exists",
				   M,
				   QMessageBox::Warning,
				   QMessageBox::Yes,
				   QMessageBox::No,
				   QMessageBox::Cancel,this);
    
    q->setButtonText( QMessageBox::Yes, "Load" );
    q->setButtonText( QMessageBox::No, "OverWrite" );
    switch(q->exec()){
    case QMessageBox::Yes:  
      Blank();
      CONTROLFILENAME=S;
      ReadChildControlFile(S);
      altered_state=false;
      Status->setText("");break;
    case QMessageBox::No:  
      CONTROLFILENAME=S;
      WriteTheOutput();
      altered_state=false;
      Status->setText("");break;
    case QMessageBox::Cancel:
      CONTROLFILENAMELine->setText(CONTROLFILENAME);  
    }
  }else{
    CONTROLFILENAME=S;
    WriteTheOutput();
    altered_state=false;
    Status->setText("");
  }
}

void Chintz_base::altered(){altered_state=true;Status->setText("Modified");}

void Chintz_base::set_CONTROLFILENAME(const QString& S){
CONTROLFILENAME=S;}

#include "generated_slotcode"
#include "generated_child_outputcode"
#include "generated_blankcode"
#include "ReadChildControlFile.cpp"
