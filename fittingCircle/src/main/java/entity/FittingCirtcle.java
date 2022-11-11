package entity;

import Jama.Matrix;

public class FittingCirtcle {
    public Matrix xc;
    public Matrix yc;
    public Matrix lastvector;

    public double getR() {
        return r;
    }

    public void setR(double r) {
        this.r = r;
    }

    public double r;
    public Matrix getXc() {
        return xc;
    }

    public void setXc(Matrix xc) {
        this.xc = xc;
    }

    public Matrix getYc() {
        return yc;
    }

    public void setYc(Matrix yc) {
        this.yc = yc;
    }

    public Matrix getLastvector() {
        return lastvector;
    }

    public void setLastvector(Matrix lastvector) {
        this.lastvector = lastvector;
    }
}
