/****************************************************************************
** Form interface generated from reading ui file 'des1.ui'
**
** Created: Mon Jun 4 14:35:57 2001
**      by:  The User Interface Compiler (uic)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/
#ifndef CHINTZ_BASE_H
#define CHINTZ_BASE_H

#include <qdialog.h>
#include <fstream>

class QVBoxLayout; 
class QHBoxLayout; 
class QGridLayout; 
class QButtonGroup;
class QCheckBox;
class QFrame;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QRadioButton;
class QSpinBox;
class QTabWidget;
class QToolButton;
class QWidget;
class QMenuBar;

class Chintz_base : public QDialog
{ 
    Q_OBJECT

public:
    Chintz_base( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~Chintz_base();

    QFrame* StatusFrame;
    QLabel* Status;
    QFrame* Framebutt;
    QPushButton* Newby;
    QPushButton* Save;
    QPushButton* OPeny;
    QMenuBar* menu;
    QLineEdit* CONTROLFILENAMELine;
    QString CONTROLFILENAME;
    QToolButton* Br;
    bool altered_state;
    ofstream ofile;

#include "generated_chintz_base_header_insert"

public slots:
    virtual void Close();
    virtual bool Clear();
    virtual void Blank();
    virtual void altered();
    virtual void WriteTheOutput();
    virtual void OpenChildControlFile();
    virtual void ChildFileSelect();
    virtual void ReadChildControlFile(const QString&);
    virtual void SaveAs();
    virtual void About();
    virtual void Help();
    virtual void control();
    virtual void WriteChildControlFile(ofstream&); 
    void set_CONTROLFILENAME(const QString&);
void RadioButton5OUTLETBOX_switch(int);
void set_TYP_BOUND(int);
void RadioButton8rando_switch(int);
void set_OPT_PT_PLACE(int);
void set_X_GRID_SIZE(const QString&);
void set_Y_GRID_SIZE(const QString&);
void set_GRID_SPACING(const QString&);
void set_NUM_PTS(const QString&);
void set_OUTLET_X_COORD(const QString&);
void set_OUTLET_Y_COORD(const QString&);
void set_SLOPED_SURF(bool);
void set_UPPER_BOUND_Z(const QString&);
void set_MEAN_ELEV(const QString&);
void RadioButton9Meshset_switch(int);
void RadioButton01TINfiles_switch(int);
void RadioButton11POfile_switch(int);
void RadioButton21ARCfile_switch(int);
void set_OPTREADINPUT(int);
void INPUTDATAFILELineFileSelect();
void set_INPUTDATAFILE(const QString&);
void set_INPUTTIME(const QString&);
void OUTFILENAMELineFileSelect();
void set_OUTFILENAME(const QString&);
void POINTFILENAMELineFileSelect();
void set_POINTFILENAME(const QString&);
void ARCGRIDFILENAMELineFileSelect();
void set_ARCGRIDFILENAME(const QString&);
void set_RUNTIME(const QString&);
void set_OPINTRVL(const QString&);
void OPTSINVARChkBxSinus_switch(int);
void set_OPTSINVAR(bool);
void OPTVARChkBxstorm_switch(int);
void set_OPTVAR(bool);
void set_SEED(const QString&);
void set_PMEAN(const QString&);
void set_STDUR(const QString&);
void set_ISTDUR(const QString&);
void set_PERIOD(const QString&);
void set_MAXPMEAN(const QString&);
void set_MAXSTDURNM(const QString&);
void set_MAXISTDURMN(const QString&);
void set_LAKEFILL(bool);
void set_INFILTRATION(const QString&);
void set_SOIL_STORE(const QString&);
void set_TRANSMISSIVITY(const QString&);
void set_KINWAVE_HQEXP(const QString&);
void set_FLOWVELOCITY(const QString&);
void set_HYDROSHAPEFAC(const QString&);
void OPTSINVARINFILTChkBxsinv_switch(int);
void set_OPTSINVARINFILT(bool);
void set_PERIOD_INFILT(const QString&);
void set_MAXICMEAN(const QString&);
void RadioButton31sinvar_switch(int);
void RadioButton31icap_switch(int);
void RadioButton41trans_switch(int);
void RadioButton51sws_switch(int);
void RadioButton51icap_switch(int);
void RadioButton61kin_switch(int);
void RadioButton71hydro_switch(int);
void set_FLOWGEN(int);
void OPTINLETChkBxInflow_switch(int);
void set_OPTINLET(bool);
void set_INLET_X(const QString&);
void set_INLET_Y(const QString&);
void set_INDRAREA(const QString&);
void set_INSEDLOAD1(const QString&);
void set_INSEDLOAD2(const QString&);
void set_INSEDLOAD3(const QString&);
void set_MAXREGDEPTH(const QString&);
void set_OPTDETACHLIM(bool);
void set_KF(const QString&);
void set_MF(const QString&);
void set_NF(const QString&);
void set_TAUCD(const QString&);
void set_PB(const QString&);
void set_OPTDIFFDEP(bool);
void set_KD(const QString&);
void set_NUMGRNSIZE(const QString&);
void set_GRAINDIAM1(const QString&);
void set_GRAINDIAM2(const QString&);
void set_GRAINDIAM3(const QString&);
void set_REGPROPORTION1(const QString&);
void set_REGPROPORTION2(const QString&);
void set_REGPROPORTION3(const QString&);
void set_BRPROPORTION1(const QString&);
void set_BRPROPORTION2(const QString&);
void set_BRPROPORTION3(const QString&);
void set_OPTREADLAYER(bool);
void set_KR(const QString&);
void set_REGINIT(const QString&);
void set_KB(const QString&);
void set_MB(const QString&);
void set_NB(const QString&);
void set_BEDROCKDEPTH(const QString&);
void set_FAULTPOS(const QString&);
void set_UPDUR(const QString&);
void set_UPRATE(const QString&);
void RadioButton81Upsy_switch(int);
void set_UPTYPE(int);
void set_BANKFULLEVENT(const QString&);
void set_HYDR_WID_COEFF_DS(const QString&);
void set_HYDR_WID_EXP_DS(const QString&);
void set_HYDR_WID_EXP_STN(const QString&);
void set_HYDR_DEP_COEFF_DS(const QString&);
void set_HYDR_DEP_EXP_DS(const QString&);
void set_HYDR_DEP_EXP_STN(const QString&);
void set_HYDR_ROUGH_COEFF_DS(const QString&);
void set_HYDR_ROUGH_EXP_DS(const QString&);
void set_HYDR_ROUGH_EXP_STN(const QString&);
void set_BANK_ROUGH_COEFF(const QString&);
void set_BANK_ROUGH_EXP(const QString&);
void OPTMNDRChkBxmndr_switch(int);
void set_OPTMNDR(bool);
void set_CRITICAL_FLOW(const QString&);
void set_DEF_CHAN_DISCR(const QString&);
void set_FRAC_WID_MOVE(const QString&);
void set_FRAC_WID_ADD(const QString&);
void set_VEG_ERODY(const QString&);
void set_LATADJUST(const QString&);
void set_BNKHTDEP(bool);
void set_OPTINTERPLAYER(bool);
void set_FP_MU(const QString&);
void set_FP_LAMBDA(const QString&);
void set_PF_DRAREAMIN(const QString&);
void set_FP_BANKFULLEVENT(const QString&);
void OPTFLOODPLAINChkBxovdep_switch(int);
void set_OPTFLOODPLAIN(bool);
void set_LOESS_DEP_RATE(const QString&);
void OPTLOESSDEPChkBxeodep_switch(int);
void set_OPTLOESSDEP(bool);
void set_OPTEXPOSURETIME(bool);
void set_BiLine_a1var(const QString&);
void set_BiLine_b1var(const QString&);
void set_BiLine_a2var(const QString&);
void set_BiLine_b2var(const QString&);
void set_BiLine_a3var(const QString&);
void set_BiLine_b3var(const QString&);
void set_BiLine_a4var(const QString&);
void set_BiLine_b4var(const QString&);
void set_BiLine_a5var(const QString&);
void set_BiLine_b5var(const QString&);
void set_BiLine_a6var(const QString&);
void set_BiLine_b6var(const QString&);
void set_BiLine_a7var(const QString&);
void set_BiLine_b7var(const QString&);
void set_BiLine_a8var(const QString&);
void set_BiLine_b8var(const QString&);
void set_BiLine_a9var(const QString&);
void set_BiLine_b9var(const QString&);
void set_BiLine_a01var(const QString&);
void set_BiLine_b01var(const QString&);
void set_BiLine_a11var(const QString&);
void set_BiLine_b11var(const QString&);
void set_BiLine_a21var(const QString&);
void set_BiLine_b21var(const QString&);

protected:
    bool event( QEvent* );
    void closeEvent( QCloseEvent* );
};

#endif // CHINTZ_BASE_H
