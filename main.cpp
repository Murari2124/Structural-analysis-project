#include<iostream>
#include<conio.h>
#include<math.h>
#include<fstream>

using namespace std;

//Function to compute each member stiffness matrix
void member_stiff_calc(float x1,float y1, float x2, float y2,int n);

//Function to assemble all member stiffness matrices
int assemble(int no,int no_mem);

//Function to multiply row and column of size n
double multiply(float m[100][100],float c[100],int n);

//Function to colve for d and f
void d_solve(int n,int);

//Function to find each member force
double member_force(int n);

//Global array of all member stiffness matrix accessible to all functions
float m_g_stiff[100][6][6];
float g_stiff[100][100];
double D[100],F[100],Q[100];


int main()
{

    float coord[50][2];
    float connection[100];
    int nnode,nmember,nDOF[50][2];

    cout<<"Welcome to the truss structure analyzer using the Stiffness method"<<endl;
    cout<<"Enter the number of nodes in the structure: ";
    cin>>nnode;
    cout<<"Enter the number of members in the structure: ";
    cin>>nmember;

    cout<<"Enter the coordinates of the node in the order you would number it.\n Also enter the DOF number associated with that node\n";


    //Loop to input the coordinates
    for(int i=0;i<nnode;i++)
    {
        cout<<"Node "<<i+1<<endl;
        cout<<"X: ";
        cin>>coord[i][0];

        cout<<"Y: ";
        cin>>coord[i][1];

        cout<<"DOF X direction: ";
        cin>>nDOF[i][0];

        cout<<"DOF Y direction: ";
        cin>>nDOF[i][1];

        cout<<endl<<endl;
    }

    //Loop to input the member conenctivity and values of A,E,L for each member
    cout<<"Now specify the connectivity between node\nPlease follow the same node numbering as before.\n Enter the node number and the required values\n";
    for(int i=0;i<nmember;i++)
    {
        cout<<"MEMBER "<<i+1;
        int near,far;
        cout<<"Enter near node number: ";
        cin>>near;
        cout<<"Enter far node number: ";
        cin>>far;

        float XN,XF,YN,YF;
        XN=coord[near-1][0];
        XF=coord[far-1][0];
        YN=coord[near-1][1];
        YF=coord[far-1][1];

        //assigning correct DOFs for ease of assembly later
        m_g_stiff[i][0][1]=m_g_stiff[i][1][0]=nDOF[near-1][0];
        m_g_stiff[i][0][2]=m_g_stiff[i][2][0]=nDOF[near-1][1];
        m_g_stiff[i][0][3]=m_g_stiff[i][3][0]=nDOF[far-1][0];
        m_g_stiff[i][0][4]=m_g_stiff[i][4][0]=nDOF[far-1][1];

        //Function to compute the member stiffness matrix in gloabal form
        member_stiff_calc(XN,YN,XF,YF,i);

        cout<<"Member stiffness matrix is: \n";
        for(int num1=0;num1<5;num1++)
        {
            for(int num2=0;num2<5;num2++)
            {
                cout<<m_g_stiff[i][num1][num2]<<"     ";
            }
            cout<<endl;
        }

    }

    //Assembling the matrix
    assemble(nnode,nmember);

    cout<<"Assembled Matrix is: ";
    for(int i=0;i<nnode*2;i++)
    {
        for(int j=0;j<nnode*2;j++)
            cout<<g_stiff[i][j]<<"   ";
        cout<<endl;
    }



    //Enter constraints and conditons
    int n_uk,n_k;
    cout<<"Number of known forces/Unknown Displacements: ";
    cin>>n_uk;
    n_k=(nnode*2)-n_uk;

    cout<<"Enter all the known forces at nodes(In order of DOF defined earlier\n";

    int num;
    for(num=1;num<=n_uk;num++)
    {
        cout<<"F"<<num<<"=";
        cin>>F[num-1];
    }

    cout<<"Enter all the known displacements/constraints/support settlements at the node\n";

    for(;num<=nnode*2;num++)
    {
        cout<<"D"<<num<<"=";
        cin>>D[num-1];
    }

    d_solve(n_uk,nnode*2);

    ofstream myfile;
    myfile.open("truss_result.csv");

    cout<<"\n\n\nALL RESULTS WILL BE SAVED IN TRUSS_RESULT.CSV(Excel File) BY THIS PROGRAM AUTOMATICALLY\n\n\n";


    myfile<<"Nodal Force "<<","<<"in N"<<endl;
    cout<<"All nodal forces are: ";
    for(int i=0;i<nnode*2;i++)
    {
        cout<<F[i]<<endl;
        myfile<<"F"<<i+1<<","<<F[i]<<endl;
    }

    myfile<<"Nodal Displacement "<<","<<"in m"<<endl;
    cout<<"All nodal displacements are: ";
    for(int i=0;i<nnode*2;i++)
    {   cout<<D[i]<<endl;
        myfile<<"D"<<i+1<<","<<D[i]<<endl;
    }

    myfile<<"Member forces "<<","<<"in N"<<endl;
    myfile<<"Positive values indicate tension and negative compression"<<endl;

    cout<<"Now member forces are: \n";
    for(int i=0;i<nmember;i++)
    {
            Q[i]=member_force(i);
            cout<<"Q"<<i+1<<Q[i];
            cout<<endl;
            myfile<<"Q"<<i+1<<","<<Q[i]<<endl;
    }


    myfile.close();

    return 0;
}

void member_stiff_calc(float x1,float y1, float x2, float y2,int n)
{
        float A,E;
        double L;
        char res;

        fflush(stdin);

            L=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
            cout<<"L="<<L<<endl;

        cout<<"Enter area of cross section (m2): ";
        cin>>A;
        cout<<"Enter E value in GPA: ";
        cin>>E;
        E=E*pow(10,9);

        float lx,ly,temp;
        lx=(x2-x1)/L;
        ly=(y2-y1)/L;
        temp=A*E/L;

        cout<<"lx= "<<lx<<endl<<"ly= "<<ly<<endl<<"AE/L= "<<temp<<endl;

        m_g_stiff[n][1][1]=m_g_stiff[n][3][3]=lx*lx;
        m_g_stiff[n][2][2]=m_g_stiff[n][4][4]=ly*ly;
        m_g_stiff[n][1][2]=m_g_stiff[n][2][1]=m_g_stiff[n][3][4]=m_g_stiff[n][4][3]=lx*ly;
        m_g_stiff[n][1][3]=m_g_stiff[n][3][1]=-1*m_g_stiff[n][1][1];
        m_g_stiff[n][2][4]=m_g_stiff[n][4][2]=-1*m_g_stiff[n][2][2];
        m_g_stiff[n][1][4]=m_g_stiff[n][2][3]=m_g_stiff[n][3][2]=m_g_stiff[n][4][1]=-1*m_g_stiff[n][1][2];

        //Storing lx*AE/l and ly*AE/l to compute member forces
        m_g_stiff[n][1][5]=lx*temp;
        m_g_stiff[n][2][5]=ly*temp;

        for(int i=1;i<=4;i++)
        {
            for(int j=1;j<=4;j++)
            {
                m_g_stiff[n][i][j]=m_g_stiff[n][i][j]*temp;
            }
        }

}

int assemble(int no,int no_mem)
{
    int ndof=no*2;

    for(int i=0;i<ndof;i++)
    {
        for(int j=0;j<ndof;j++)
        {
            g_stiff[i][j]=0;

            //Loop to search all member matrices
            for(int k=0;k<no_mem;k++)
            {
                for(int l=1;l<=4;l++)
                {
                    if(m_g_stiff[k][l][0]==(i+1))
                    {
                        for(int m=1;m<=4;m++)
                        {
                            if(m_g_stiff[k][0][m]==(j+1))
                            {   g_stiff[i][j]=g_stiff[i][j]+m_g_stiff[k][l][m];
                                break;
                            }
                        }
                    }
                }

            }
        }
    }
    return 0;
}

/*double multiply(float m[100][100],float c[100],int n)
{
 	double sum=0;
	for(int i=0;i<n;i++)
	{
		cout<<m[n][i]<<"    "<<c[i]<<endl;
		sum=sum+m[n][i]*c[i];
	}
	cout<<sum;
	return sum;
}
*/
void d_solve(int n,int n_tot)
{
    //Inversion of matrix
    float matrix[100][100],inv[100][100],t,a;
    int i, j, k;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            matrix[i][j]=g_stiff[i][j];
        }
    }


    for(i = 0; i < n; i++)
    {
        for(j = n; j < 2*n; j++)
        {
            if(i==(j-n))
                matrix[i][j] = 1.0;
            else
                matrix[i][j] = 0.0;
        }
    }

    for(i=0;i<n;i++)
    {
        t=matrix[i][i];
        for(j=i;j<2*n;j++)
            matrix[i][j]=matrix[i][j]/t;
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
                t=matrix[j][i];
                for(k=0;k<2*n;k++)
                    matrix[j][k]=matrix[j][k]-t*matrix[i][k];
            }
        }
    }

    printf("The inverse matrix is: \n");
    for(i = 0; i < n; i++){
        for(j = n,k=0; j < 2*n; j++,k++)
        {
            cout<<matrix[i][j];
            inv[i][k]=matrix[i][j];
            printf("\t");
        }
        printf("\n");
    }


    for(i=0;i<n; i++)
    {
        D[i]=0.0;
        for(int j=0;j<n;j++)
        {
            //cout<<inv[i][j]<<"     "<<F[j]<<endl;
            D[i]=D[i]+(inv[i][j]*F[j]);
        }

    }

    for(;i<n_tot;i++)
    {
        F[i]=0.0;
        for(int j=0;j<n;j++)
        {
            //cout<<inv[i][j]<<"     "<<F[j]<<endl;
            F[i]=F[i]+(g_stiff[i][j]*F[j]);
        }

    }
}

double member_force(int n)
{

    float v1,v2;
    double sum;
    v1= m_g_stiff[n][1][5];
    v2= m_g_stiff[n][2][5];

    double temp1[4],temp2[4];

    temp1[0]=-1*v1;
    temp1[1]=-1*v2;
    temp1[2]=v1;
    temp1[3]=v2;

    int index[4];
    index[0]=m_g_stiff[n][0][1];
    index[1]=m_g_stiff[n][0][2];
    index[2]=m_g_stiff[n][0][3];
    index[3]=m_g_stiff[n][0][4];

    temp2[0]=D[index[0]-1];
    temp2[1]=D[index[1]-1];
    temp2[2]=D[index[2]-1];
    temp2[3]=D[index[3]-1];

    sum=0.0;
    for(int i=0;i<4;i++)
        sum=sum+(temp1[i]*temp2[i]);

    return sum;
}
