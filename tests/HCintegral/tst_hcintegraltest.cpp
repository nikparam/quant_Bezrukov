#include <QString>
#include <QtTest>

class HCintegralTest : public QObject
{
    Q_OBJECT

public:
    HCintegralTest();

private Q_SLOTS:
    void testCase1();
};

HCintegralTest::HCintegralTest()
{
}

void HCintegralTest::testCase1()
{
    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(HCintegralTest)

#include "tst_hcintegraltest.moc"
