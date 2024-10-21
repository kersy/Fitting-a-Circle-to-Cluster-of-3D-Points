package entity;

public class Ellipse {
    //椭圆方程参数构成：x^2 + A*x*y + B*y^2 + C*x + D*y + E = 0;
    double A ;//椭圆方程参数A
    double B ;//椭圆方程参数B
    double C ;//椭圆方程参数C
    double D ;//椭圆方程参数D
    double E ;//椭圆方程参数E
    double Ea ;//椭圆长轴
    double Eb ;//椭圆短轴
    double ova;//椭圆度，单位为千分之数值

    public double getEa() {
        return Ea;
    }

    public void setEa(double a) {
        this.Ea = a;
    }
    public double getEb() {
        return Eb;
    }

    public void setEb(double b) {
        this.Eb = b;
    }
    public double getOva() {
        return ova;
    }

    public void setOva(double ova) {
        this.ova = ova;
    }

    public Ellipse() {

    }


    public double getA() {
        return A;
    }

    public void setA(double a) {
        A = a;
    }

    public double getB() {
        return B;
    }

    public void setB(double b) {
        B = b;
    }

    public double getC() {
        return C;
    }

    public void setC(double c) {
        C = c;
    }

    public double getD() {
        return D;
    }

    public void setD(double d) {
        D = d;
    }

    public double getE() {
        return E;
    }

    public void setE(double e) {
        E = e;
    }

}
