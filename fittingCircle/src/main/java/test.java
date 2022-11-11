import Jama.Matrix;
import Jama.SingularValueDecomposition;
import entity.FittingCirtcle;
public class test {
    public static double norm(Matrix x)
    {
        double [][] t = (x.arrayTimes(x)).getArray();
        double num=0;
        for(int i =0;i< x.getRowDimension();i++)
        {
            num+=t[i][0];
            num+=t[i][1];
            num+=t[i][2];
        }
        return Math.sqrt(num);
    }
    public static Matrix cross(Matrix n0,Matrix n1,int sig)//向量积运算
    {
        if(sig==0)
        {
            n0=n0.times(-1);
        }
        else if (sig==1){
            n0=n0.times(1);
        }
        double[][] n0t = n0.getArray();
        double[][] n1t = n1.getArray();
        double[][] mat = new double[2][3];
        mat[0][0]=n0t[0][0];
        mat[0][1]=n0t[0][1];
        mat[0][2]=n0t[0][2];
        mat[1][0]=n1t[0][0];
        mat[1][1]=n1t[0][1];
        mat[1][2]=n1t[0][2];
        Matrix detMatrix = new Matrix(mat);
        Matrix ni= detMatrix.getMatrix(0,1,1,2);
        int[] row={0,1};
        int[] column = {0,2};
        Matrix nj= detMatrix.getMatrix(row,column);
        Matrix nk= detMatrix.getMatrix(0,1,0,1);
        double nin=ni.det();
        double nkn=nk.det();
        double njn=-1*nj.det();
        double[][] K ={{nin,njn,nkn}};
        Matrix k =new Matrix(K);
        return k;
    }
    public static FittingCirtcle fit_circle_2d(Matrix P_xy)//二维平面内拟合圆
    {
        Matrix one = new Matrix(P_xy.getRowDimension(),1,1);
        Matrix x =P_xy.getMatrix(0,P_xy.getRowDimension()-1,0,0);
        Matrix y =P_xy.getMatrix(0,P_xy.getRowDimension()-1,1,1);
        P_xy.setMatrix(0,P_xy.getRowDimension()-1,2,2,one);
        Matrix A=P_xy.transpose();
        Matrix b= x.arrayTimes(x).plus(y.arrayTimes(y));
        A=A.transpose();
        Matrix c = A.solve(b);
        double[][] ct = c.getArray();
        double xc = ct[0][0]/2;
        double yc = ct[1][0]/2;
        double[][] lastvector={{xc},{yc},{0}};
        double r = Math.sqrt(ct[2][0]+xc*xc+yc*yc);
        FittingCirtcle fittingCirtcle = new FittingCirtcle();
        Matrix xct = new Matrix(1,1,xc);
        Matrix yct = new Matrix(1,1,yc);
        Matrix lstv = new Matrix(lastvector);
        fittingCirtcle.setLastvector(lstv);
        fittingCirtcle.setXc(xct);
        fittingCirtcle.setYc(yct);
        fittingCirtcle.setR(r);
        return fittingCirtcle;
    }
    public static  Matrix rodrigues_rot (Matrix P,Matrix n0,Matrix n1,int sig)//罗德里格斯旋转
    {
        n0=n0.times(1/norm(n0));
        n1=n1.times(1/norm(n1));
        Matrix k = cross(n0,n1,sig);
        if(k.norm1()!=0)
        {

            k=k.times(1/norm(k));
        }
        Matrix theta=n0.times(n1.transpose());
        double[][] thetatemp=theta.getArray();
        double thetanum =Math.acos(thetatemp[0][0]);
        double[][] Pront = new double[P.getRowDimension()][3];
        for(int i=0;i<P.getRowDimension();i++)
        {
            Matrix pi = P.getMatrix(i,i,0,2);
            Matrix t1 = pi.times(Math.cos(thetanum));
            Matrix t2 = cross(k,pi,1).times(Math.sin(thetanum));
            Matrix t3 = k.times(k.times(pi.transpose()).getArray()[0][0]).times(1-Math.cos(thetanum));
            Matrix temp = t1.plus(t2).plus(t3);
            Pront[i][0]=temp.getArray()[0][0];
            Pront[i][1]=temp.getArray()[0][1];
            Pront[i][2]=temp.getArray()[0][2];
        }
        Matrix P_rot=new Matrix(Pront);
        return P_rot;
    }
    public static void main(String[] args) {
        double[][] array = {
//                {2,0,0},
//                {0,2,0},
//                {0,-2,0},
                {197297.4617, 216355.2787, 448.6076},
                {197293.7126, 216349.8803, 451.8861},
                {197296.2875, 216353.5825, 451.7551},
                {197294.8214, 216351.4768, 452.4308},
                {197294.9019, 216351.6648, 443.4885},
                {197296.7171, 216354.3012, 444.8446},
                {197296.92, 216354.601, 445.2778},
                {197292.6999, 216348.4749, 445.701},
                {197293.65, 216349.8288, 444.0722},
                {197292.3469, 216347.9819, 447.7407},
                {197297.3921, 216355.3633, 448.1248}


        };
        Matrix P=new Matrix(array);
        double[][]line1= P.getArray();
        double row1mean = 0;
        double row2mean = 0;
        double row3mean = 0;
        for(int i=0;i<line1.length;i++)
        {
            row1mean+=line1[i][0];

            row2mean+=line1[i][1];

            row3mean+=line1[i][2];
        }
        double[][] p_meantemp = new double[line1.length][3];
        row1mean=row1mean/line1.length;
        row2mean=row2mean/line1.length;
        row3mean=row3mean/line1.length;
        for(int i=0;i<line1.length;i++)
        {
            p_meantemp[i][0]=row1mean;
            p_meantemp[i][1]=row2mean;
            p_meantemp[i][2]=row3mean;

        }

        Matrix P_MEAN=new Matrix(p_meantemp);
        double[][] tempmean = {{row1mean},{row2mean},{row3mean}};
        Matrix P_mean = new Matrix(tempmean);
        Matrix P_center = P.minus(P_MEAN);
        SingularValueDecomposition singularValueDecomposition =  P_center.svd();
        Matrix U = singularValueDecomposition.getU();
        Matrix s = singularValueDecomposition.getS();
        Matrix V = singularValueDecomposition.getV();
        V=V.inverse();
        V=V.getMatrix(2,2,0,2);
        Matrix normal = V;
        normal=normal.times(-1);
        double[][] n1temp={{0,0,1}};
        Matrix n1 = new Matrix(n1temp);
        Matrix P_xy = rodrigues_rot(P_center,normal,n1,1);
        FittingCirtcle circle = fit_circle_2d(P_xy);
        Matrix laststep = circle.getLastvector().transpose();
        Matrix C = rodrigues_rot(laststep,n1,normal,1).plus(P_mean.transpose());
        System.out.println("fitting circlecenter:");
        C.print(1,7);
        System.out.println("r="+circle.r);;
        System.out.println("法向量n：");
        normal.print(1,6);
        //最终结果为Matrix 调用getarray方法得到二维数组取得其中的值
        //奇异值分解后的正负号差别影响最后N法向量的方向，不影响描述具体平面
        //所用范数为第二范数，即根号下矩阵元素平方和
    }

}
