#include "chintz_base.h"
#include <qapplication.h>
#include <qmainwindow.h>

int main(int argc, char** argv)
{
    QApplication app(argc,argv);
    Chintz_base c;

    if ( QApplication::desktop()->width() > c.width() + 20
	 && QApplication::desktop()->height() > c.height() +20 )
	c.show();
    else
	c.showMaximized();

    c.show();

    qApp->setMainWidget((QWidget*) &c);

    return app.exec();
}
