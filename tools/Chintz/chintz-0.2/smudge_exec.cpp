#include <fstream>
#include <string>

ifstream ifl("smudge.ui");
ofstream ofl("generated_chintz_base_constructor_insert");
ofstream hed("generated_chintz_base_header_insert");
ofstream slo("generated_slot_header_insert");
ofstream sloco("generated_slotcode");
ofstream output("generated_child_outputcode");
ofstream blank("generated_blankcode");
ofstream reader("generated_child_readcode");
ofstream bireader("generated_extra_readcode");

string Str;

void write(string type,string parent,string enabled,string name,string title,string bold,
	   string x0,string y0,string x1,string y1){
  ofl<<name<<"->setEnabled( "<<enabled<<" );"<<endl;
  ofl<<name<<"->setGeometry( QRect( "<<x0<<", "<<y0<<", "<<x1<<", "<<y1<<" ) );"<<endl; 
  ofl<<"QFont "<<name<<"_font="<<name<<"->font();"<<endl;
  ofl<<name<<"_font.setBold( "<<bold<<" );"<<endl;
  ofl<<name<<"->setFont("<<name<<"_font);"<<endl;

  if (type=="Label" || type=="RadioButton" || type=="CheckBox")
    ofl<<name<<"->setText( tr( "<<'"'<<title<<'"'<<"  ) );"<<endl;

  if (type=="Box" || type=="Buttons" && title!="")ofl<<name<<"->setTitle( tr( "<<'"'<<title<<'"'<<"  ) );"<<endl;
  if (type=="Frame"){
    if (title!=""){
      ofl<<"QLabel* "<<name<<"Label=new QLabel("<<name<<','<<'"'<<title<<'"'<<"  );"<<endl;
      ofl<<name<<"Label->setGeometry(20,10,100,20);"<<endl;
      ofl<<name<<"Label->setText( tr( "<<'"'<<title<<'"'<<"  ) );"<<endl;
      ofl<<"QFont "<<name<<"Lfont="<<name<<"Label->font();"<<endl;
      ofl<<name<<"Lfont.setItalic( true );"<<endl;
      ofl<<name<<"Lfont.setBold( false );"<<endl;
      ofl<<name<<"Lfont.setPointSize( 12 );"<<endl;
      ofl<<name<<"->setFont("<<name<<"Lfont);"<<endl;
      ofl<<"QPalette "<<name<<"pal;"<<endl;
      ofl<<"QColorGroup "<<name<<"cg;"<<endl;
      ofl<<name<<"cg.setColor( QColorGroup::Foreground, QColor( 0, 0, 255) );"<<endl;
      ofl<<name<<"cg.setColor( QColorGroup::Background, QColor( 192, 192,192) );"<<endl;
      ofl<<name<<"pal.setActive("<<name<<"cg);"<<endl;
      ofl<<name<<"pal.setInactive("<<name<<"cg);"<<endl;
      ofl<<name<<"Label->setPalette("<<name<<"pal);"<<endl;
    }
  }


    
}
void Chomp(string& Title,bool& found){
  unsigned s;
  if ((s=Str.find("<Title>"))==0){
    found=true;Str.erase(0,7);
    Title="";
    if (Str.find("<")==string::npos)Title=Str+" ";
    while ((s=Str.find("<"))==string::npos){ifl>>Str;if (Str.find("<")==string::npos)Title=Title+Str+" ";}
    if ((s=Str.find("<"))!=0){Title=Title+Str.substr(0,s);Str.erase(0,s);}
  }
}

void Slurp(string& b){
  unsigned c;
  if (Str.length()==0)ifl>>Str;
  if ((c=Str.find("<"))!=string::npos){
    b=Str.substr(0,c);
    Str.erase(0,c);
  }else{
    b=Str;Str="";
  }
}
class Line{
  static int num;
public:
  static void make(string parent){
    string Name="TextLabel";
    string Lr="Line";
    string Br;
    string Title="Title";
    string x0="10";
    string y0="60";
    string x1="16";
    string y1="20";
    string xx0="30";
    string yy0="60";
    string xx1="40";
    string yy1="22";
    string Labelbold="FALSE";
    string Linebold="FALSE";
    string Lineenable="TRUE";
    
    string state;
    string var;
    bool found;

    unsigned s;
    
    string id="0123456789",h;
    int r,t;
    num++;
    r=num;
    while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
    Name=Name+h;
    Title=Title+h;
    Lr=Lr+h;
    var=Lr+"var";

    if ((s=Str.find("id="))==0){int w=Str.find(">");Lr=Str.substr(3,w-3);Str.erase(0,Lr.length()+4);var=Lr;Lr=Lr+"Line";}

    cout<<Lr<<endl; 
 
    if (Str.length()==0)ifl>>Str;
    while (!ifl.eof()){
     found=false; 
     Chomp(Title,found);
     if (found) state="Title";
      if (state=="Title"){
	if (Str.find("</Title>")==0){found=true;Str.erase(0,8);state="";}
	if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Labelbold);}

	if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
      }else{
	if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(xx0);ifl>>yy0>>xx1;Slurp(yy1);}
	if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Linebold);}
	if (Str.find("<Enabled>")==0){found=true;Str.erase(0,9);Slurp(Lineenable);}
      }
      if(Str.find("<BrowseButton>")==0){
	found=true;Str.erase(0,15);
	Br=Lr+"Tool";
	hed<<"QToolButton* "<<Br<<";"<<endl;
	ofl<<Br<<" = new QToolButton( "<<parent<<", "<<'"'<<Br<<'"'<<" );"<<endl;
	ofl<<Br<<"->setGeometry( QRect( "<<atoi(xx0.c_str())+atoi(xx1.c_str())<<", "<<atoi(yy0.c_str())<<", 20, 20 ) );"<<endl; 
	ofl<<Br<<"->setText( tr( "<<'"'<<"..."<<'"'<<"  ) );"<<endl;
	ofl<<"QFont "<<Br<<"_font="<<Br<<"->font();"<<endl;
	ofl<<Br<<"_font.setBold( false );"<<endl;
	ofl<<Br<<"->setFont("<<Br<<"_font);"<<endl;
	ofl<<"QToolTip::add(  "<<Br<<", tr( "<<'"'<<"browse"<<'"'<<" ) );"<<endl;
	slo<<"void "<<Lr<<"FileSelect();"<<endl;
	sloco<<"void Chintz_base::"<<Lr<<"FileSelect(){"<<endl;
	sloco<<"QString S=QFileDialog::getSaveFileName(QString::null,"<<'"'<<"*"<<'"'<<",this);"<<endl;
	sloco<<"if (S.length()!=0){set_"<<var<<"(S);"<<Lr<<"->setText(S);}}"<<endl;
        ofl<<"connect("<<Br<<", SIGNAL( clicked() ), this, SLOT( "<<Lr<<"FileSelect() ) );"<<endl;

      }
      if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
      if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
      if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
      if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
      if (Str.find("</BrowseButton>")==0){found=true;Str.erase(0,16);}
      if (!found && Str.find("</Line>")!=0){
	cout<<"Smudge::Tag Error in <Line>! Offending string is: "<<Str<<endl;
	exit(1);
      }
      if (Str.length()==0)ifl>>Str;
      if (Str.find("</Line>")==0){Str.erase(0,7);break;}
    }
    hed<<"QLabel* "<<Name<<";"<<endl;
    hed<<"QLineEdit* "<<Lr<<";"<<endl;

    ofl<<Name<<" = new QLabel( "<<parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
    ofl<<Lr<<" = new QLineEdit( "<<parent<<", "<<'"'<<Lr<<'"'<<" );"<<endl;

    write("Label",parent,Lineenable,Name,Title,Labelbold,x0,y0,x1,y1);
    write("Line",parent,Lineenable,Lr,"",Linebold,xx0,yy0,xx1,yy1);
    
    slo<<"void set_"<<var<<"(const QString&);"<<endl;
 
    hed<<"QString "<<var<<";"<<endl;
    sloco<<"void Chintz_base::set_"<<var<<"(const QString& S){"<<endl;
    sloco<<var<<"=S;}"<<endl;

    output<<"if ("<<var<<".length() != 0){"<<endl;
    output<<"  ofile<<"<<'"'<<var<<":"<<'"'<<"<<endl;"<<endl;
    output<<"  ofile<<"<<var<<"<<endl;"<<endl;
    output<<"}"<<endl;

    blank<<Lr<<"->setText("<<'"'<<'"'<<");"<<endl;
    blank<<"set_"<<var<<"("<<'"'<<'"'<<");"<<endl;

    reader<<"if(val.find("<<'"'<<var<<":"<<'"'<<")!=string::npos){"<<endl; 
    reader<<"detected=true;"<<endl;
    reader<<Lr<<"->setText(qty);"<<endl;
    reader<<"set_"<<var<<"(qty);}"<<endl;

    ofl<<"connect("<<Lr<<", SIGNAL( textChanged(const QString&) ), this, SLOT( set_"<<var<<"(const QString&) ) );"<<endl;
     ofl<<"connect("<<Lr<<", SIGNAL( textChanged(const QString&) ), this, SLOT( altered() ) );"<<endl;
   }
};
int Line::num=0;

class BiLine{
  static int num;
public:
  static void make(string parent){
    string Name="BiLine_a";
    string Lr="BiLine_b";
    string Br;
    string Title="Title";
    string x0="10";
    string y0="60";
    string x1="150";
    string y1="20";
    string xx0="30";
    string yy0="60";
    string xx1="40";
    string yy1="22";
    string Labelbold="FALSE";
    string Linebold="FALSE";
    string Lineenable="TRUE";
    
    string state;
    string var;
    bool found;
    unsigned s;
    
    string id="0123456789",h;
    int r,t;
    num++;
    r=num;
    while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
    Name=Name+h;
    Title=Title+h;
    Lr=Lr+h;

    cout<<Lr<<endl; 
 
    if (Str.length()==0)ifl>>Str;
    while (!ifl.eof()){
      found=false;
      if (Str.find("<Origin>")==0){found=true;Str.erase(0,8);Slurp(x0);Slurp(y0);}
      if (Str.find("</Origin>")==0){found=true;Str.erase(0,9);}
      if (!found && Str.find("</BiLine>")!=0){
	cout<<"Smudge::Tag Error in <BiLine>! Offending string is: "<<Str<<endl;
	exit(1);
      }
    
      if (Str.length()==0)ifl>>Str;
      if (Str.find("</BiLine>")==0){found=true;Str.erase(0,9);break;}
    }
    hed<<"QLineEdit* "<<Name<<";"<<endl;
    hed<<"QLineEdit* "<<Lr<<";"<<endl;

    ofl<<Name<<" = new QLineEdit( "<<parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
    ofl<<Lr<<" = new QLineEdit( "<<parent<<", "<<'"'<<Lr<<'"'<<" );"<<endl;

    ofl<<"QLabel* "<<Name<<"Label=new QLabel("<<parent<<","<<'"'<<Name<<"Label"<<'"'<<");"<<endl;
    ofl<<""<<Name<<"Label->setText(tr("<<'"'<<"Variable Name"<<'"'<<"));"<<endl;
    ofl<<""<<Name<<"Label->setGeometry(QRect("<<x0<<','<<y0<<','<<x1<<','<<y1<<"));"<<endl;

    ofl<<"QLabel* "<<Lr<<"Label=new QLabel("<<parent<<","<<'"'<<Lr<<"Label"<<'"'<<");"<<endl;
    ofl<<""<<Lr<<"Label->setText(tr("<<'"'<<"Value"<<'"'<<"));"<<endl;
    ofl<<""<<Lr<<"Label->setGeometry(QRect("<<x0<<','<<atoi(y0.c_str())+2*atoi(y1.c_str())<<','<<x1<<','<<y1<<"));"<<endl;

    ofl<<Name<<"->setGeometry("<<x0<<','<<atoi(y0.c_str())+atoi(y1.c_str())<<','<<x1<<','<<y1<<");"<<endl;
    ofl<<Lr<<"->setGeometry("<<x0<<','<<atoi(y0.c_str())+3*atoi(y1.c_str())<<','<<x1<<','<<y1<<");"<<endl;
    
    slo<<"void set_"<<Name<<"var(const QString&);"<<endl;
    slo<<"void set_"<<Lr<<"var(const QString&);"<<endl;
 
    hed<<"QString "<<Name<<"var;"<<endl;
    hed<<"QString "<<Lr<<"var;"<<endl;

    sloco<<"void Chintz_base::set_"<<Name<<"var(const QString& S){"<<endl;
    sloco<<Name<<"var"<<"=S;}"<<endl;

    sloco<<"void Chintz_base::set_"<<Lr<<"var(const QString& S){"<<endl;
    sloco<<Lr<<"var"<<"=S;}"<<endl;

    output<<"if ("<<Name<<"var.length() != 0 && "<<Lr<<"var.length()!=0){"<<endl;
    output<<"  ofile<<"<<Name<<"var<<"<<'"'<<":"<<'"'<<"<<endl;"<<endl;
    output<<"  ofile<<"<<Lr<<"var<<endl;"<<endl;
    output<<"}"<<endl;

    blank<<Name<<"->setText("<<'"'<<'"'<<");"<<endl;
    blank<<"set_"<<Name<<"var"<<"("<<'"'<<'"'<<");"<<endl;

    blank<<Lr<<"->setText("<<'"'<<'"'<<");"<<endl;
    blank<<"set_"<<Lr<<"var"<<"("<<'"'<<'"'<<");"<<endl;

    bireader<<"if(num=="<<num<<"){"<<endl;
    bireader<<"int pos=val.find("<<'"'<<":"<<'"'<<");"<<endl;
    bireader<<"string f=val.substr(0,pos);"<<endl; 
    bireader<<Name<<"->setText(f.c_str());"<<endl;
    bireader<<"set_"<<Name<<"var(f.c_str());"<<endl;
    bireader<<Lr<<"->setText(qty);"<<endl;
    bireader<<"set_"<<Lr<<"var(qty);}"<<endl;

  ofl<<"connect("<<Name<<", SIGNAL( textChanged(const QString&) ), this, SLOT( set_"<<Name<<"var"<<"(const QString&) ) );"<<endl;
     ofl<<"connect("<<Name<<", SIGNAL( textChanged(const QString&) ), this, SLOT( altered() ) );"<<endl;

  ofl<<"connect("<<Lr<<", SIGNAL( textChanged(const QString&) ), this, SLOT( set_"<<Lr<<"var"<<"(const QString&) ) );"<<endl;
     ofl<<"connect("<<Lr<<", SIGNAL( textChanged(const QString&) ), this, SLOT( altered() ) );"<<endl;
   }
};
int BiLine::num=0;

class Box{
  static int num;
public:
  static void make(string);
};
class Buttons{
  static int num;
public:
  static void make(string);
};
class RadioButton{
  static int num;
public:
  static void make(string);
};
class CheckBox{
  static int num;
public:
  static void CheckBox::make(string);
};
int CheckBox::num=0;

class Frame{
  static int num;
public:
  static void make(string Parent){
    bool found;
    unsigned s;
    string Name="Frame",Title="",Bold="TRUE",Enabled="TRUE",Style="StyledPanel",Shadow="Raised";
    string x0="10",y0="20",x1="375",y1="260";
    string id="0123456789",h;
    int r,t;
    num++;
    r=num;
    while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
    Name=Name+h;
    //Title=Title+h;
    if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);}
    cout<<Name<<endl; 

    hed<<"QFrame* "<<Name<<";"<<endl;
    ofl<<Name<<" = new QFrame( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;

    if (Str.length()==0)ifl>>Str;
    while (!ifl.eof()){
      found=false;
      if ((s=Str.find("<Frame"))==0){found=true;Str.erase(0,7);Frame::make(Name);}
      if (Str.find("<Enabled>")==0){found=true;Str.erase(0,9);Slurp(Enabled);}
      if ((s=Str.find("<Box"))==0){found=true;Str.erase(0,5);Box::make(Name);}
      if ((s=Str.find("<CheckBox"))==0){found=true;Str.erase(0,10);CheckBox::make(Name);}
      if ((s=Str.find("<Buttons"))==0){found=true;Str.erase(0,9);Buttons::make(Name);}
      if ((s=Str.find("<Line"))==0){found=true;Str.erase(0,6);Line::make(Name);}
      if ((s=Str.find("<BiLine"))==0){found=true;Str.erase(0,8);BiLine::make(Name);}
      if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
      if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}

      Chomp(Title,found);

      if ((s=Str.find("<Style>"))==0){
	found=true;Str.erase(0,7);Slurp(Style);
      }
      if ((s=Str.find("<Shadow>"))==0){
	found=true;Str.erase(0,8);Slurp(Shadow);
      }
      if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
      if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
      if (Str.find("</Style>")==0){found=true;Str.erase(0,8);}
      if (Str.find("</Shadow>")==0){found=true;Str.erase(0,9);}
      if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
      if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
      if (!found && Str.find("</Frame>")!=0){
	cout<<"Smudge::Tag Error in <Frame>! Offending string is: "<<Str<<endl;
	exit(1);
      }
      if (Str.length()==0)ifl>>Str;
      if (Str.find("</Frame>")==0){found=true;Str.erase(0,8);break;}
    }
    write("Frame",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
    ofl<<Name<<"->setFrameShape( QFrame::"<<Style<<" );"<<endl;
    ofl<<Name<<"->setFrameShadow( QFrame::"<<Shadow<<" );"<<endl;
  }
};
int Frame::num=0;

void Box::make(string Parent){
  unsigned s;
  bool found;
  string Name="Box",Title="",Bold="TRUE",Enabled="TRUE";
  string x0="10",y0="20",x1="100",y1="100";
  string id="0123456789",h;
  int r,t;
  num++;
  r=num;
  while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
  Name=Name+h;
  
  if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);}
    cout<<Name<<endl; 

  hed<<"QGroupBox* "<<Name<<";"<<endl;
  ofl<<Name<<" = new QGroupBox( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
  
  if (Str.length()==0)ifl>>Str;
  while (!ifl.eof()){
    found=false;
    if ((s=Str.find("<Frame"))==0){found=true;Str.erase(0,7);Frame::make(Name);}
    if ((s=Str.find("<CheckBox"))==0){found=true;Str.erase(0,10);CheckBox::make(Name);}
    if ((s=Str.find("<Buttons"))==0){found=true;Str.erase(0,9);Buttons::make(Name);}
    if ((s=Str.find("<Box"))==0){found=true;Str.erase(0,5);Box::make(Name);}
    if (Str.find("<Enabled>")==0){found=true;Str.erase(0,9);Slurp(Enabled);}
    if ((s=Str.find("<Line"))==0){found=true;Str.erase(0,6);Line::make(Name);}
    if ((s=Str.find("<BiLine"))==0){found=true;Str.erase(0,8);BiLine::make(Name);}
    if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}
    if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
 
    Chomp(Title,found);
 
    if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
    if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
    if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
    if (!found && Str.find("</Box>")!=0){
      cout<<"Smudge::Tag Error in <Box>! Offending string is: "<<Str<<endl;
      exit(1);
    }
    if (Str.length()==0)ifl>>Str;
    if (Str.find("</Box>")==0){Str.erase(0,6);break;}
  }
  write("Box",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
}
int Box::num=0;

void Buttons::make(string Parent){
  unsigned s;
  bool found;
  string Name="Buttons",Title="Buttons",Bold="TRUE",Enabled="TRUE",var;
  string x0="10",y0="20",x1="100",y1="100";
  string id="0123456789",h;
  int r,t;
  num++;
  r=num;
  while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
  Name=Name+h;
  Title=Title+h;
  
  var=Name+"var";

  if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);var=Name;Name=Name+"BGroup";}

  cout<<Name<<endl;

  hed<<"QButtonGroup* "<<Name<<";"<<endl;
  ofl<<Name<<" = new QButtonGroup( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
  ofl<<Name<<"->setExclusive(true);"<<endl;
  
  if (Str.length()==0)ifl>>Str;
  while (!ifl.eof()){
    found=false;
    if ((s=Str.find("<RadioButton>"))==0){found=true;Str.erase(0,13);RadioButton::make(Name);}
    if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
    if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}

    Chomp(Title,found);

    if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
    if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
    if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
    if (!found && Str.find("</Buttons>")!=0){
      cout<<"Smudge::Tag Error in <Buttons>! Offending string is: "<<Str<<endl;
      exit(1);
    }
    if (Str.length()==0)ifl>>Str;
    if (Str.find("</Buttons>")==0){Str.erase(0,10);break;}
  }
  write("Buttons",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
  ofl<<Name<<"->setFrameShadow( QButtonGroup::Raised );"<<endl;

  hed<<"QString "<<var<<";"<<endl;
  
  slo<<"void set_"<<var<<"(int);"<<endl;
  
  sloco<<"void Chintz_base::set_"<<var<<"(int i){"<<endl;
  sloco<<" "<<var<<".setNum(i);"<<endl;
  sloco<<"}"<<endl;
  
  ofl<<"connect("<<Name<<", SIGNAL( clicked(int) ), this, SLOT( set_"<<var<<"(int) ) );"<<endl;
  ofl<<"connect("<<Name<<", SIGNAL( clicked(int) ), this, SLOT( altered() ) );"<<endl;

  blank<<Name<<"->setButton(0);"<<endl;
  blank<<"set_"<<var<<"(0);"<<endl;

  reader<<"if(val.find("<<'"'<<var<<":"<<'"'<<")!=string::npos){"<<endl; 
  reader<<"detected=true;"<<endl;
  reader<<Name<<"->setButton(qty.toInt());"<<endl;
  reader<<"set_"<<var<<"(qty.toInt());}"<<endl;

  output<<"if ("<<var<<".length() != 0){"<<endl;
  output<<"  ofile<<"<<'"'<<var<<":"<<'"'<<"<<endl;"<<endl;
  output<<"  ofile<<"<<var<<"<<endl;"<<endl;
  output<<"}"<<endl;
}
int Buttons::num=0;

void RadioButton::make(string Parent){
  unsigned s;
  bool found;
  string Name="RadioButton",Title="RadioButton",Bold="TRUE",Enabled="TRUE",Checked="FALSE";
  string x0="10",y0="20",x1="100",y1="20";
  string widge,poo;
  string id="0123456789",h;
  int r,t;
  bool swoo=false;

  num++;
  r=num;
  while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
  Name=Name+h;
  Title=Title+h;
  
  if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);}

    cout<<Name<<endl; 

  hed<<"QRadioButton* "<<Name<<";"<<endl;
  ofl<<Name<<" = new QRadioButton( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
  
  if (Str.length()==0)ifl>>Str;
  while (!ifl.eof()){
    found=false;
    if ((s=Str.find("<RadioButton>"))==0){found=true;Str.erase(0,13);RadioButton::make(Name);}
    if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
    if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}
    if (Str.find("<Switch>")==0){
      found=true;Str.erase(0,8);
      Slurp(widge);Slurp(poo);
      slo<<"void "<<Name<<widge<<"_switch(int);"<<endl;

      sloco<<"void Chintz_base::"<<Name<<widge<<"_switch(int i){"<<endl;
      sloco<<"if("<<widge<<"==0)return;"<<endl;
      sloco<<"if (i !=0) {"<<widge <<"->setEnabled("<<poo<<");}"<<endl;
      sloco<<" else {"<<widge <<"->setEnabled(!"<<poo<<");}}"<<endl;

      ofl<<"connect("<<Name<<", SIGNAL( stateChanged(int) ), this, SLOT( "<<Name<<widge<<"_switch(int) ) );"<<endl;
    }

    Chomp(Title,found);

    if ((s=Str.find("<Checked>"))==0){
      found=true;Str.erase(0,9);Slurp(Checked);
    }
    if (Str.find("</Switch>")==0){found=true;Str.erase(0,9);}
    if (Str.find("</Checked>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
    if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
    if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
    if (!found && Str.find("</RadioButton>")!=0){
      cout<<"Smudge::Tag Error in <RadioButton>! Offending string is: "<<Str<<endl;
      exit(1);
    }
    if (Str.length()==0)ifl>>Str;
    if (Str.find("</RadioButton>")==0){Str.erase(0,14);break;}
  }
  write("RadioButton",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
  ofl<<Name<<"->setChecked( "<<Checked<<" );"<<endl;
}
int RadioButton::num=0;

  
void CheckBox::make(string Parent){
  unsigned s;
  bool found;
  string Name="CheckBox",Title="RadioButton",Bold="TRUE",Enabled="TRUE",Checked="FALSE",var;
  string x0="10",y0="20",x1="100",y1="100";
  string id="0123456789",h;
  string widge,poo;
  int r,t;
  num++;
  r=num;
  while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
  Name=Name+h;
  Title=Title+h;

  var=Name+"var";
   if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);var=Name;Name=Name+"ChkBx";}

  hed<<"QCheckBox* "<<Name<<";"<<endl;
  ofl<<Name<<" = new QCheckBox( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
  
  if (Str.length()==0)ifl>>Str;
  while (!ifl.eof()){
    found=false;
    if ((s=Str.find("<CheckBox"))==0){found=true;Str.erase(0,10);CheckBox::make(Name);}
    if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
    if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}
    if (Str.find("<Switch>")==0){
      found=true;Str.erase(0,8);
      Slurp(widge);Slurp(poo);
      slo<<"void "<<Name<<widge<<"_switch(int);"<<endl;

      sloco<<"void Chintz_base::"<<Name<<widge<<"_switch(int i){"<<endl;
      sloco<<"if("<<widge<<"==0)return;"<<endl;
      sloco<<"if (i !=0) {"<<widge <<"->setEnabled("<<poo<<");}"<<endl;
      sloco<<" else {"<<widge <<"->setEnabled(!"<<poo<<");}}"<<endl;

      ofl<<"connect("<<Name<<", SIGNAL( stateChanged(int) ), this, SLOT( "<<Name<<widge<<"_switch(int) ) );"<<endl;
    }

    Chomp(Title,found);

    if ((s=Str.find("<Checked>"))==0){
      found=true;Str.erase(0,9);Slurp(Checked);
    }
    if (Str.find("</Switch>")==0){found=true;Str.erase(0,9);}
    if (Str.find("</Checked>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
    if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
    if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
    if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
    if (!found && Str.find("</CheckBox>")!=0){
      cout<<"Smudge::Tag Error in <CheckBox>! Offending string is: "<<Str<<endl;
      exit(1);
    }
    if (Str.length()==0)ifl>>Str;
    if (Str.find("</CheckBox>")==0){Str.erase(0,11);break;}
  }
  write("CheckBox",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
  ofl<<Name<<"->setChecked( "<<Checked<<" );"<<endl;

  hed<<"QString "<<var<<";"<<endl;
  
  slo<<"void set_"<<var<<"(bool);"<<endl;
  

  blank<<Name<<"->setChecked(0);"<<endl;
  blank<<"set_"<<var<<"(0);"<<endl;

  reader<<"if(val.find("<<'"'<<var<<":"<<'"'<<")!=string::npos){"<<endl; 
  reader<<"detected=true;"<<endl;
  reader<<Name<<"->setChecked(qty.toInt());"<<endl;
  reader<<"set_"<<var<<"(qty.toInt());}"<<endl;

  sloco<<"void Chintz_base::set_"<<var<<"(bool i){"<<endl;
  sloco<<" "<<var<<".setNum(i);"<<endl;
  sloco<<"}"<<endl;
  
  ofl<<"connect("<<Name<<", SIGNAL( toggled(bool) ), this, SLOT( set_"<<var<<"(bool) ) );"<<endl;
  ofl<<"connect("<<Name<<", SIGNAL( clicked() ), this, SLOT( altered() ) );"<<endl;

  output<<"if ("<<var<<".length() != 0){"<<endl;
  output<<"  ofile<<"<<'"'<<var<<":"<<'"'<<"<<endl;"<<endl;
  output<<"  ofile<<"<<var<<"<<endl;"<<endl;
  output<<"}"<<endl;
  }


class Page{
  static int num;
public:
  static void make(string Parent){
    unsigned s;
    bool found;
    string Name="Page",Title="Page",Bold="TRUE",Enabled="TRUE";
    string x0="10",y0="50",x1="750",y1="520";
    string id="0123456789",h;
    int r,t;
    num++;
    r=num;
    while (r>0){t=r;r=r/10;t=t-r*10;h=h+id[t];}
    Name=Name+h;
    Title=Title+h;

    if ((s=Str.find("id="))==0){int w=Str.find(">");Name=Str.substr(3,w-3);Str.erase(0,Name.length()+4);}
    //        cout<<Name<<endl; 

    hed<<"QWidget* "<<Name<<";"<<endl;
    ofl<<Name<<" = new QWidget( "<<Parent<<", "<<'"'<<Name<<'"'<<" );"<<endl;
    if (Str.length()==0)ifl>>Str;

    while (!ifl.eof()){
      found=false;
      Chomp(Title,found);
      if ((s=Str.find("<Frame"))==0){found=true;Str.erase(0,7);Frame::make(Name);}
      if ((s=Str.find("<Box"))==0){found=true;Str.erase(0,5);Box::make(Name);}
      if ((s=Str.find("<Line"))==0){found=true;Str.erase(0,6);Line::make(Name);}
      if ((s=Str.find("<BiLine"))==0){found=true;Str.erase(0,8);BiLine::make(Name);}
      if (Str.find("<Bold>")==0){found=true;Str.erase(0,6);Slurp(Bold);}
      if (Str.find("<Geometry>")==0){found=true;Str.erase(0,10);Slurp(x0);ifl>>y0>>x1;Slurp(y1);}
      if (Str.find("</Bold>")==0){found=true;Str.erase(0,7);}
      if (Str.find("</Enabled>")==0){found=true;Str.erase(0,10);}
      if (Str.find("</Title>")==0){found=true;Str.erase(0,8);}
      if (Str.find("</Geometry>")==0){found=true;Str.erase(0,11);}
      if (!found && Str.find("</Page>")!=0){
	cout<<"Smudge::Tag Error in <Page>! Offending string is: "<<Str<<endl;
	exit(1);
      }
      if (Str.length()==0)ifl>>Str;
      if (Str.find("</Page>")==0){Str.erase(0,7);break;}
    }
    write("Page",Parent,Enabled,Name,Title,Bold,x0,y0,x1,y1);
    ofl<<Parent<<"->insertTab( "<<Name<<", tr( "<<'"'<<Title<<'"'<<" ) );"<<endl;

  }
};
int Page::num=0;

int main(){
  string Name="Tab";
  unsigned s;
  

  output<<"void Chintz_base::WriteChildControlFile(ofstream& ofile){"<<endl;
  output<<"ofile<<"<<'"'<<"#Automatically generated CHILD input file Chintz v0.2"<<'"'<<"<<endl;"<<endl;

  blank<<"void Chintz_base::Blank(){"<<endl;

  hed<<"QTabWidget* "<<Name<<";"<<endl;
  ofl<<Name<<" = new QTabWidget( this, "<<'"'<<Name<<'"'<<" );"<<endl;
  ofl<<Name<<"->setGeometry( QRect( 10, 70, 770, 520 ) );"<<endl; 
  


  ifl>>Str;
  while (!ifl.eof()){
    if ((s=Str.find("<Page"))==0){Str.erase(s,6);Page::make(Name);}
    if (Str.length()==0)ifl>>Str;
  }


  output<<"}"<<endl;
  blank<<"}"<<endl;
}
