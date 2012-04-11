#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QtGui>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->openButton, SIGNAL(clicked()), this, SLOT(OpenImage()));
    connect(ui->saveButton, SIGNAL(clicked()), this, SLOT(SaveImage()));
    connect(ui->harrisButton, SIGNAL(clicked()), this, SLOT(HarrisCornerImage()));
    connect(ui->matchButton, SIGNAL(clicked()), this, SLOT(MatchImages()));
    connect(ui->RANSACButton, SIGNAL(clicked()), this, SLOT(RANSAC()));
    connect(ui->stitchButton, SIGNAL(clicked()), this, SLOT(StitchImages()));

    ui->harrisSpinBox->setValue(2.0);
    ui->harrisThresSpinBox->setValue(50.0);
    ui->RANSACThresSpinBox->setValue(5.0);
    ui->iterationsBox->setValue(200);

    ui->tabWidget->setCurrentIndex(0);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::DrawDisplayImage()
{
    ui->img1Display->setPixmap(QPixmap::fromImage(m_InImage1Display));
    ui->img2Display->setPixmap(QPixmap::fromImage(m_InImage2Display));
    ui->imgSDisplay->setPixmap(QPixmap::fromImage(m_StitchedImage));
}



void MainWindow::OpenImage()
{
    const QString title;

    QString fileName = QFileDialog::getOpenFileName(this, title);

    if(ui->tabWidget->currentIndex() == 0)
    {
        if (!fileName.isEmpty())
            m_InImage1.load(fileName);

        m_InImage1Display = m_InImage1.copy();
    }

    if(ui->tabWidget->currentIndex() == 1)
    {
        if (!fileName.isEmpty())
            m_InImage2.load(fileName);

        m_InImage2Display = m_InImage2.copy();
    }

    DrawDisplayImage();

}

void MainWindow::SaveImage()
{
    const QString title;

    QString fileName = QFileDialog::getSaveFileName(this, title);

    if(ui->tabWidget->currentIndex() ==  0)
    {
        if (!fileName.isEmpty())
            m_InImage1Display.save(fileName);
    }

    if(ui->tabWidget->currentIndex() ==  1)
    {
        if (!fileName.isEmpty())
            m_InImage2Display.save(fileName);
    }

    if(ui->tabWidget->currentIndex() ==  2)
    {
        if (!fileName.isEmpty())
            m_StitchedImage.save(fileName);
    }
}

void MainWindow::HarrisCornerImage()
{
    double sigma = ui->harrisSpinBox->value();
    double thres = ui->harrisThresSpinBox->value();

    HarrisCornerDetector(m_InImage1, sigma, thres, &m_IntPts1, m_NumIntPts1, m_InImage1Display);
    HarrisCornerDetector(m_InImage2, sigma, thres, &m_IntPts2, m_NumIntPts2, m_InImage2Display);

    DrawDisplayImage();
}

void MainWindow::MatchImages()
{
    m_InImage1Display = m_InImage1.copy();
    m_InImage2Display = m_InImage2.copy();

    MatchInterestPoints(m_InImage1, m_IntPts1, m_NumIntPts1,
               m_InImage2, m_IntPts2, m_NumIntPts2,
               &m_Matches, m_NumMatches, m_InImage1Display, m_InImage2Display);

    DrawDisplayImage();
}

void MainWindow::RANSAC()
{
    int numIterations = ui->iterationsBox->value();
    double inlierThreshold = ui->RANSACThresSpinBox->value();

    m_InImage1Display = m_InImage1.copy();
    m_InImage2Display = m_InImage2.copy();

    RANSAC(m_Matches, m_NumMatches, numIterations, inlierThreshold,
           m_Hom, m_HomInv, m_InImage1Display, m_InImage2Display);

    DrawDisplayImage();
}

void MainWindow::StitchImages()
{
    Stitch(m_InImage1, m_InImage2, m_Hom, m_HomInv, m_StitchedImage);

    ui->tabWidget->setCurrentIndex(2);
    DrawDisplayImage();
}





