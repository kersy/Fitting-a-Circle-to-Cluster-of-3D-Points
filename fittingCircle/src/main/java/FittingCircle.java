import Jama.Matrix;
import Jama.SingularValueDecomposition;
import entity.Ellipse;
import entity.FittingCirtcle;
public class FittingCircle {
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

    public static Ellipse EllipseFitting(double[][] Ecount)//椭圆拟合算法
    {
        double x2y2=0;
        double xy3=0;
        double x2y=0;
        double xy2=0;
        double xy=0;

        double y4=0;
        double y3=0;
        double y2=0;

        double x2=0;
        double x=0;

        double y=0;
        double N=Ecount.length;

        double x3y=0;
        double x3=0;


        for(int i =0 ;i<Ecount.length;i++)
        {
            double xt =Ecount[i][1];
            double yt =Ecount[i][2];
            x2y2+=xt*xt*yt*yt;
            xy3+=xt*yt*yt*yt;
            x2y+=xt*xt*yt;
            xy2+=xt*yt*yt;
            xy+=xt*yt;

            y4+=yt*yt*yt*yt;
            y3+=yt*yt*yt;
            y2+=yt*yt;

            x2+=xt*xt;
            x+=xt;

            y+=yt;


            x3y+=xt*xt*xt*yt;
            x3+=xt*xt*xt;
        }
        double[][] M1t={
                {x2y2,xy3,x2y,xy2,xy},
                {xy3,y4,xy2,y3,y2},
                {x2y,xy2,x2,xy,x},
                {xy2,y3,xy,y2,y},
                {xy,y2,x,y,N}
        };
        double[][] M2t={
                {-x3y},
                {-x2y2},
                {-x3},
                {-x2y},
                {-x2}
        };
        Matrix M1 = new Matrix(M1t);
        Matrix M2 = new Matrix(M2t);
    try{
        Matrix M0 = M1.inverse().times(M2);
        System.out.println(M1.rank());
        double[][]M0t = M0.getArray();
        double A = M0t[0][0];
        double B = M0t[1][0];
        double C = M0t[2][0];
        double D = M0t[3][0];
        double E = M0t[4][0];

        double Xp = (A*D-2*B*C)/(A*A-4*B);
        double Yp = (A*C-2*D)/(A*A-4*B);
        double Xc = -Xp;
        double Yc = -Yp;
        double theta_r = 0.5*Math.atan(A/(B-1));
        double theta_offset = -theta_r;
        double a = Math.sqrt((Xp*Xp+A*Xp*Yp+B*Yp*Yp-E)/(Math.pow(Math.cos(theta_r),2)-A*Math.sin(theta_r)*Math.cos(theta_r)+B*Math.pow(Math.sin(theta_r),2)));
        double b = Math.sqrt((Xp*Xp+A*Xp*Yp+B*Yp*Yp-E)/(Math.pow(Math.sin(theta_r),2)+A*Math.sin(theta_r)*Math.cos(theta_r)+B*Math.pow(Math.cos(theta_r),2)));
        Ellipse ellipse = new Ellipse();
        ellipse.setA(A);
        ellipse.setB(B);
        ellipse.setC(C);
        ellipse.setD(D);
        ellipse.setE(E);//写入椭圆参数
        ellipse.setEa(Math.max(a,b)*2);
        ellipse.setEb(Math.min(a,b)*2);
        ellipse.setOva((ellipse.getEa()-ellipse.getEb())/((ellipse.getEa()+ellipse.getEb()))*2000);//写入椭圆度，已化为千分之计量
        return ellipse;
    }catch(Exception e)
        {
            System.out.println("奇异矩阵拟合失败");
        }
    return null;
    }
    public static void main(String[] args) {
        double[][] array = {
//                {0,1,0},
//                {0,2,0},
//                {0,-2,0},
//                {197297.4617, 216355.2787, 448.6076},
//                {197293.7126, 216349.8803, 451.8861},
//                {197296.2875, 216353.5825, 451.7551},
//                {197294.8214, 216351.4768, 452.4308},
//                {197294.9019, 216351.6648, 443.4885},
//                {197296.7171, 216354.3012, 444.8446},
//                {197296.92, 216354.601, 445.2778},
//                {197292.6999, 216348.4749, 445.701},
//                {197293.65, 216349.8288, 444.0722},
//                {197292.3469, 216347.9819, 447.7407},
//                {197297.3921, 216355.3633, 448.1248}
//                {	0	,	47836.7812	,	292.4664	},
//                {	0	,	47836.4604	,	293.8825	},
//                {	0	,	47836.3276	,	295.2579	},
//                {	0	,	47836.3381	,	296.7034	},
//                {	0	,	47836.4695	,	298.1027	},
//                {	0	,	47836.6907	,	299.217	},
//                {	0	,	47837.0016	,	300.0415	},
//                {	0	,	47837.8443	,	300.0772	},
//                {	0	,	47838.0825	,	299.4701	},
//                {	0	,	47838.4458	,	295.4761	},
//                {	0	,	47838.3051	,	293.9825	},
//                {	0	,	47838.0693	,	292.8761	},
//                {	0	,	47837.7679	,	292.1229	},
//                {	0	,	47837.2675	,	291.7947	},

                {65544.4491,-3.0636,17.7779},
                {65544.4423,-2.5206,18.6527},
                {65544.455,-1.8904,19.3795},
                {65544.474,-0.6602,20.3182},
                {65544.4713,0.6029,20.8798},
                {65544.4792,2.4475,21.1716},
                {65544.4693,3.6362,21.0669},
                {65544.465,4.545,20.8252},
                {65546.4608,4.6946,20.7399},
                {65546.466,3.6442,21.0373},


        };
        Matrix P=new Matrix(array);
        if (P.getRowDimension() == 3) {
            Matrix n1;
            Matrix n2;
            Matrix n3;
            System.out.println(P.det());
            if (Math.abs((int) P.det()) == 0) {
                n1 = P.getMatrix(0, 0, 0, 2);
                n2 = P.getMatrix(1, 1, 0, 2);
                n3 = P.getMatrix(2, 2, 0, 2);
                int n1t = (int) Math.abs(n1.times(n2.transpose()).times(1 / (norm(n1) * norm(n2))).get(0, 0));
                int n2t = (int) Math.abs(n2.times(n3.transpose()).times(1 / (norm(n2) * norm(n3))).get(0, 0));
                int n3t = (int) Math.abs(n1.times(n3.transpose()).times(1 / (norm(n1) * norm(n3))).get(0, 0));

                if (Math.acos(n1t) == 0 && Math.acos(n2t) == 0 && Math.acos(n3t) == 0) {

                    return;
                }

            }
        }
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
        //椭圆拟合
        double[][] n2temp={{1,0,0}};//椭圆拟合向量 投影至yoz平面
        Matrix n2 = new Matrix(n2temp);
        Matrix P_yz = rodrigues_rot(P,normal,n2,1);//椭圆拟合向量
        System.out.println("rodYZ");
        P_yz.print(1,2);
        Ellipse ellipse = EllipseFitting(P_yz.getArray());//椭圆拟合完成
        //计算椭长短轴、圆度并输出
        System.out.println(ellipse.getEa());//输出长轴
        System.out.println(ellipse.getEb());//输出短轴
        System.out.println(ellipse.getOva());//输出椭圆度
    }

}
