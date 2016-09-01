

//read a child input file 
void Chintz_base::ReadChildControlFile(const QString& inputfile){
  ifstream infile;
  infile.open(inputfile);
  string q,val;
  QString qty;
  char c;
  bool detected;
  int num=0;
  while (!infile.eof())
    {
      infile>>val;
      if (val.find("Comments")==0)break;
      detected=false;
      if(!infile.good())break;
      infile.get(c);
      while (c!='\n')infile.get(c);
      if(val.find("#")==0)continue;
      q="#";
      while(q.find("#")==0)infile>>q;
      infile.get(c);
      while (c!='\n')infile.get(c);
      qty=q.c_str();
      if (qty.isNull()) continue;
      if (qty.isEmpty()) continue;
      #include "generated_child_readcode"
      if (!detected){
       num++;
       #include "generated_extra_readcode"
      }
    }
}

