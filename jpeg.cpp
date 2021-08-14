#include <iostream>
#include <stdio.h>
#include <cmath>
#include <malloc.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <bitset>

const double PI = 3.141592653589793;
const double SQH = 0.707106781186547;  /* square root of 2 */

using namespace std;

string code="";
string press="";

static void fft1(float *data,int nn,int isign);

inline int bin2int(string bin)
{
   int sum = 0;
   for( int i=0; i<bin.length(); i++ )
       sum = (sum<<1) | (bin[i]&1);
   return sum;
}

void dct1(float *x,int n)
{
    int i,ii,nn,mm;
    float tc,ts,sqn,temp2;
    double temp1;
    float *v;

    nn = n >> 1;
    mm = n << 1;
    sqn = (float)sqrt((double)n);

    v = (float *) calloc (mm,sizeof(float));
    if (v == NULL)
    {
        printf("allocation failure\n");
        exit(1);
    }

    for (i=0; i<nn; i++)
    {
        ii = i << 1;
        v[ii] = x[ii];
        v[ii+1] = 0.0;
    }
    for (i=nn; i<n; i++)
    {
        ii = i << 1;
        v[ii] = x[mm-ii-1];
        v[ii+1] = 0.0;
    }

    fft1(v-1,n,1);

    temp2 = SQH/sqn;
    x[0] = v[0]/sqn;
    for (i=1; i<=nn; i++)
    {
        ii = i << 1;
        temp1 = (double)(PI*i/mm);
        tc = (float) cos(temp1);
        ts = (float) sin(temp1);
        x[i] = 2.0*(tc*v[ii] + ts*v[ii+1])*temp2;
        x[n-i] = 2.0*(ts*v[ii] - tc*v[ii+1])*temp2;
    }

    free(v);
}

void idct1(float *x,int n)
{
    int i,ii,mm,nn;
    float *v;
    float temp2,tc,ts,sqn;
    double temp1;

    nn = n >> 1;
    mm = n << 1;
    sqn = (float)sqrt((double)n);

    v = (float *) calloc (mm,sizeof(float));
    if (v == NULL)
    {
        printf("allocation failure\n");
        exit(1);
    }

    temp2 = sqn/SQH;
    v[0] = x[0]*sqn;
    v[1] = 0.0;
    for (i=1; i<n; i++)
    {
        ii = i << 1;
        temp1 = (double)(PI*i/mm);
        tc = (float)cos(temp1);
        ts = (float)sin(temp1);
        v[ii] = 0.5*(tc*x[i] + ts*x[n-i])*temp2;
        v[ii+1] = 0.5*(ts*x[i] - tc*x[n-i])*temp2;
    }

    fft1(v-1,n,-1);

    for (i=0; i<nn; i++)
    {
        ii = i << 1;
        x[ii] = v[ii];
    }
    for (i=nn; i<n; i++)
    {
        ii = i << 1;
        x[mm-ii-1] = v[ii];
    }
    free(v);
}

static void fft1(float *data,int nn,int isign)
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;
    n = nn << 1;
    j = 1;
    for (i=1; i<n; i+=2)
    {
        if (j>i)
        {
            swap(data[j],data[i]);
            swap(data[j+1],data[i+1]);
        }
        m = n >> 1;
        while (m>=2 && j>m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n>mmax)
    {
        istep = 2*mmax;
        theta = 6.28318530717959/(-isign*mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m<mmax; m+=2)
        {
            for (i=m; i<=n; i+=istep)
            {
                j = i+mmax;
                tempr = wr*data[j]-wi*data[j+1];
                tempi = wr*data[j+1]+wi*data[j];
                data[j] = data[i]-tempr;
                data[j+1] = data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp=wr)*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }

    if (isign == -1)
    {
        for (i=1; i<=n; ++i)
            data[i] = data[i]/nn;
    }
}

void dct2(float **x,int n)
{
    int i,j;
    float *y;

    y = (float *) calloc (n,sizeof(float));
    if (y == NULL)
    {
        printf("allocation failure\n");
        exit(1);
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            y[j] = x[j][i];
        dct1(y,n);
        for (j=0; j<n; j++)
            x[j][i] = y[j];
    }   /* end of loop i */

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            y[j] = x[i][j];
        dct1(y,n);
        for (j=0; j<n; j++)
            x[i][j] = y[j];
    }   /* end of loop i */

    free(y);
}

void idct2(float **x,int n)
{
    int i,j;
    float *y;

    y = (float *) calloc (n,sizeof(float));
    if (y == NULL)
    {
        printf("allocation failure\n");
        exit(1);
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            y[j] = x[j][i];
        idct1(y,n);
        for (j=0; j<n; j++)
            x[j][i] = y[j];
    }   /* end of loop i */

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            y[j] = x[i][j];
        idct1(y,n);
        for (j=0; j<n; j++)
            x[i][j] = y[j];
    }   /* end of loop i */

    free(y);
}

float quant[8][8]={16, 11, 10, 16, 24, 40, 51, 61,
                   12, 12, 14, 19, 26, 58, 60, 55,
                   14, 13, 16, 24, 40, 57, 69, 56,
                   14, 17, 22, 29, 51, 87, 80, 62,
                   18, 22, 37, 56, 68, 109, 103, 77,
                   24, 36, 55, 64, 81, 104, 113, 92,
                   49, 64, 78, 87, 103, 121, 120, 101,
                   72, 92, 95, 98, 112, 100, 103, 99};

void fixQuant(int coe)
{
    if(coe>=50)
        coe=200-2*coe;
    else
        coe=5000/coe;

    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            quant[i][j]*=((float)coe/100);
            if(quant[i][j]>255)
                quant[i][j]=255;
            else if(quant[i][j]<=0)
                quant[i][j]=1;
        }
    }
}

void quantize(float **input)
{
    int turn=0;
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            turn=input[i][j]/quant[i][j]+0.5;
            input[i][j]=turn;
        }
    }
}

unsigned int dcCode[12]={0b00,0b010,0b011,0b100,0b101,0b110,0b1110,0b11110,0b111110,0b1111110,0b11111110,0b111111110};

int dcCodeLen[12]={2,3,3,3,3,3,4,5,6,7,8,9};

void dcencode(float first, float ex)
{
    int diff = 0;
    diff = first-ex;
    int check=abs(diff);
    int power=0;
    while(check!=0)
    {
        check/=2;
        power++;
    }
    string mld="";
    mld = bitset<9>(dcCode[power]).to_string();
    mld.erase(0,9-dcCodeLen[power]);
    code+=mld;

    string word1="";
    if(diff==0)
        word1+="0";
    else if(diff>0)
    {
        while(diff!=0)
        {
            if(diff%2==0)
                word1+="0";
            else
                word1+="1";
            diff/=2;
        }
    }
    else
    {
        int diff2=abs(diff);
        while(diff2!=0)
        {
            if(diff2%2==1)
                word1+="0";
            else
                word1+="1";
            diff2/=2;
        }
    }
    reverse(word1.begin(),word1.end());
    code+=word1;
}

int order[63]={1 , 8 , 16, 9 , 2 , 3 , 10, 17,
               24, 32, 25, 18, 11, 4 , 5 , 12,
               19, 26, 33, 40, 48, 41, 34, 27,
               20, 13, 6 , 7 , 14, 21, 28, 35,
               42, 49, 56, 57, 50, 43, 36, 29,
               22, 15, 23, 30, 37, 44, 51, 58,
               59, 52, 45, 38, 31, 39, 46, 53,
               60, 61, 54, 47, 55, 62, 63};

unsigned int acCode[16][10]={
    0b00,0b01,0b100,0b1011,0b11010,0b1111000,0b11111000,0b1111110110,0b1111111110000010,0b1111111110000011,
    0b1100,0b11011,0b1111001,0b111110110,0b11111110110,0b1111111110000100,0b1111111110000101,0b1111111110000110,0b1111111110000111,0b1111111110001000,
    0b11100,0b11111001,0b1111110111,0b111111110100,0b1111111110001001,0b1111111110001010,0b1111111110001011,0b1111111110001100,0b1111111110001101,0b1111111110001110,
    0b111010,0b111110111,0b111111110101,0b1111111110001111,0b1111111110010000,0b1111111110010001,0b1111111110010010,0b1111111110010011,0b1111111110010100,0b1111111110010101,
    0b111011,0b1111111000,0b1111111110010110,0b1111111110010111,0b1111111110011000,0b1111111110011001,0b1111111110011010,0b1111111110011011,0b1111111110011100,0b1111111110011101,
    0b1111010,0b11111110111,0b1111111110011110,0b1111111110011111,0b1111111110100000,0b1111111110100001,0b1111111110100010,0b1111111110100011,0b1111111110100100,0b1111111110100101,
    0b1111011,0b111111110110,0b1111111110100110,0b1111111110100111,0b1111111110101000,0b1111111110101001,0b1111111110101010,0b1111111110101011,0b1111111110101100,0b1111111110101101,
    0b11111010,0b111111110111,0b1111111110101110,0b1111111110101111,0b1111111110110000,0b1111111110110001,0b1111111110110010,0b1111111110110011,0b1111111110110100,0b1111111110110101,
    0b111111000,0b111111111000000,0b1111111110110110,0b1111111110110111,0b1111111110111000,0b1111111110111001,0b1111111110111010,0b1111111110111011,0b1111111110111100,0b1111111110111101,
    0b111111001,0b1111111110111110,0b1111111110111111,0b1111111111000000,0b1111111111000001,0b1111111111000010,0b1111111111000011,0b1111111111000100,0b1111111111000101,0b1111111111000110,
    0b111111010,0b1111111111000111,0b1111111111001000,0b1111111111001001,0b1111111111001010,0b1111111111001011,0b1111111111001100,0b1111111111001101,0b1111111111001110,0b1111111111001111,
    0b1111111001,0b1111111111010000,0b1111111111010001,0b1111111111010010,0b1111111111010011,0b1111111111010100,0b1111111111010101,0b1111111111010110,0b1111111111010111,0b1111111111011000,
    0b1111111010,0b1111111111011001,0b1111111111011010,0b1111111111011011,0b1111111111011100,0b1111111111011101,0b1111111111011110,0b1111111111011111,0b1111111111100000,0b1111111111100001,
    0b11111111000,0b1111111111100010,0b1111111111100011,0b1111111111100100,0b1111111111100101,0b1111111111100110,0b1111111111100111,0b1111111111101000,0b1111111111101001,0b1111111111101010,
    0b1111111111101011,0b1111111111101100,0b1111111111101101,0b1111111111101110,0b1111111111101111,0b1111111111110000,0b1111111111110001,0b1111111111110010,0b1111111111110011,0b1111111111110100,
    0b1111111111110101,0b1111111111110110,0b1111111111110111,0b1111111111111000,0b1111111111111001,0b1111111111111010,0b1111111111111011,0b1111111111111100,0b1111111111111101,0b1111111111111110
};

int acCodeLen[16][10]={
     2, 2, 3, 4, 5, 7, 8,10,16,16,
     4, 5, 7, 9,11,16,16,16,16,16,
     5, 8,10,12,16,16,16,16,16,16,
     6, 9,12,16,16,16,16,16,16,16,
     6,10,16,16,16,16,16,16,16,16,
     7,11,16,16,16,16,16,16,16,16,
     7,12,16,16,16,16,16,16,16,16,
     8,12,16,16,16,16,16,16,16,16,
     9,15,16,16,16,16,16,16,16,16,
     9,16,16,16,16,16,16,16,16,16,
     9,16,16,16,16,16,16,16,16,16,
    10,16,16,16,16,16,16,16,16,16,
    10,16,16,16,16,16,16,16,16,16,
    11,16,16,16,16,16,16,16,16,16,
    16,16,16,16,16,16,16,16,16,16,
    16,16,16,16,16,16,16,16,16,16
};

void acencode(float **go)
{
    vector <int> rlc (0);
    int zerocount=0;
    for(int i=0;i<63;i++)
    {
        if(go[order[i]/8][order[i]%8]!=0)
        {
            rlc.push_back(zerocount);
            rlc.push_back(go[order[i]/8][order[i]%8]);
            zerocount=0;
        }
        else
        {
            zerocount++;
            if(zerocount==15)
            {
                rlc.push_back(-1);
                zerocount=0;
            }
        }
    }
    rlc.push_back(-2);

    int counting=0;
    while(counting < rlc.size())
    {
        if(rlc[counting]==-2)
            code+="1010";
        else if(rlc[counting]==-1)
        {
            bool eob=false;
            int zrl=0;
            while(rlc[counting]<0)
            {
                if(rlc[counting]==-2)
                {
                    code+="1010";
                    eob=true;
                    counting=rlc.size()+10;
                    break;
                }
                zrl++;
                counting++;
            }

            if(!eob)
            {
                for(int i=0;i<zrl;i++)
                    code+="11111111001";
                counting-=1;
            }
        }
        else
        {
            int num=abs(rlc[counting+1]);
            int power2=0;
            while(num!=0)
            {
                num/=2;
                power2++;
            }

            if(rlc[counting]==0&&power2==1)
                code+="00";
            else if(rlc[counting]==0&&power2==2)
                code+="01";
            else
            {
                string mid="";
                mid = bitset<16>(acCode[rlc[counting]][power2-1]).to_string();
                while(mid[0]=='0')
                    mid.erase(0,1);
                code+=mid;
            }

            string word2="";
            if(rlc[counting+1]>0)
            {
                int b=rlc[counting+1];
                while(b!=0)
                {
                    if(b%2==0)
                        word2+="0";
                    else
                        word2+="1";
                    b/=2;
                }
            }
            else
            {
                int c=abs(rlc[counting+1]);
                while(c!=0)
                {
                    if(c%2==1)
                        word2+="0";
                    else
                        word2+="1";
                    c/=2;
                }
            }
            reverse(word2.begin(), word2.end());
            code+=word2;
            counting++;
        }
        counting++;
    }
}

struct chunk
{
    float **block=new float*[8];
};

void dcdecode(chunk backing)
{
    string category="";
    bool fin=false;
    int next=0;
    while(!fin)
    {
        category+=press[0];
        press.erase(0, 1);
        for(int i=0;i<12;i++)
        {
            if(category.length()==dcCodeLen[i])
            {
                if(bin2int(category)==dcCode[i])
                {
                    next=i;
                    fin=true;
                    break;
                }
            }
            else if(category.length()<dcCodeLen[i])
                break;
        }
    }

    category="";
    if(next==0)
    {
        press.erase(0, 1);
        backing.block[0][0]=0;
        return;
    }
    else
        category.assign(press, 0, next);

    int dcValue=0;
    if(category[0]=='1')
    {
        for(int i=0;i<category.length();i++)
        {
            dcValue*=2;
            if(category[i]=='1')
                dcValue+=1;
        }
    }
    else
    {
        for(int i=0;i<category.length();i++)
        {
            dcValue*=2;
            if(category[i]=='0')
                dcValue+=1;
        }
        dcValue*=-1;
    }
    press.erase(0, next);
    backing.block[0][0]=dcValue;
}

void acdecode(chunk backing)
{
    string rscode="";
    bool fin=false;
    int run=0;
    int coe=0;
    int now=0;
    while(!fin)
    {
        rscode+=press[0];
        press.erase(0, 1);
        if(rscode=="1010")
        {
            fin=true;
            while(now!=63)
            {
                backing.block[order[now]/8][order[now]%8]=0;
                now++;
            }
        }
        else if(rscode=="11111111001")
        {
            for(int i=0;i<15;i++)
            {
                backing.block[order[now]/8][order[now]%8]=0;
                now++;
            }
            rscode="";
        }
        else
        {
            for(int i=0;i<16;i++)
            {
                for(int j=0;j<10;j++)
                {
                    if(rscode.length()<acCodeLen[i][j])
                        break;
                    else if(rscode.length()==acCodeLen[i][j])
                    {
                        if(bin2int(rscode)==acCode[i][j])
                        {
                            run=i;
                            for(int a=0;a<run;a++)
                            {
                                backing.block[order[now]/8][order[now]%8]=0;
                                now++;
                            }

                            coe=j+1;
                            rscode.assign(press, 0, coe);
                            int acValue=0;
                            if(rscode[0]=='1')
                            {
                                for(int a=0;a<rscode.length();a++)
                                {
                                    acValue*=2;
                                    if(rscode[a]=='1')
                                        acValue+=1;
                                }
                            }
                            else
                            {
                                for(int a=0;a<rscode.length();a++)
                                {
                                    acValue*=2;
                                    if(rscode[a]=='0')
                                        acValue+=1;
                                }
                                acValue*=-1;
                            }
                            backing.block[order[now]/8][order[now]%8]=acValue;
                            now++;
                            press.erase(0, coe);
                            rscode="";
                            i=20;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void deQuantize(float **input)
{
    for(int i=0;i<8;i++)
    {
        for(int j=0;j<8;j++)
        {
            input[i][j]*=quant[i][j];
        }
    }
}

string show(char c)
{
    bitset<8> btmp(c);
    return btmp.to_string();
}

int main ()
{
    int operation=0;
    cout << endl << "*Choose an operation:(1 to compress/2 to depress):";
    cin >> operation;

    string filename;
    cout << endl << "*Input filename:";
    cin >> filename;

    string outfilename;
    if(operation==1||operation==2)
    {
        cout << endl << "*Input output filename:";
        cin >> outfilename;
    }

    int qf=0;
    cout << endl << "*Input quality factor:";
    cin >> qf;

    fixQuant(qf);

    if(operation==1)
    {
        int height=512, width=512;
        char image[height][width];

        FILE *f0 = fopen(filename.c_str(), "rb");
        fread(image, sizeof(unsigned char), width*height, f0);
        fclose(f0);

        vector <chunk> pic;
        for(int i=0;i<64;i++)
        {
            for(int j=0;j<64;j++)
            {
                chunk temp;
                for(int a=0;a<8;a++)
                {
                    temp.block[a] = new float[8];
                    for(int b=0;b<8;b++)
                    {
                        temp.block[a][b]=(unsigned char)image[i*8+a][j*8+b];
                        temp.block[a][b]-=128;
                    }
                }
                pic.push_back(temp);
            }
        }

        for(int i=0;i<4096;i++)
        {
            dct2(pic[i].block, 8);
            quantize(pic[i].block);
        }

        dcencode(pic[0].block[0][0], 0);
        acencode(pic[0].block);
        for(int i=1;i<4096;i++)
        {
            dcencode(pic[i].block[0][0], pic[i-1].block[0][0]);
            acencode(pic[i].block);
        }

        FILE *f1 = fopen(outfilename.c_str(), "wb");
        int eight=0;
        unsigned char output;
        for(int i=0;i<code.size();i++)
        {
            output<<=1;
            if(code[i]=='1')
                output+=1;
            eight++;
            if(eight==8)
            {
                fwrite(&output, sizeof(char), 1, f1);
                eight=0;
                output=0;
            }
        }
        while(eight!=8)
        {
            output<<=1;
            eight++;
        }
        fwrite(&output, sizeof(char), 1, f1);
    }
    else if(operation==2)
    {
        FILE *f1 = fopen(filename.c_str(), "rb");
        char buff;
        while (!feof(f1))
        {
            fread(&buff, sizeof(char), 1, f1);
            press+=show(buff);
            for (int i = 0; i < 8; i++)
                buff <<= 1;
        }

        vector <chunk> reduct;
        for(int i=0;i<4096;i++)
        {
            chunk temp;
            for(int j=0;j<8;j++)
                temp.block[j] = new float[8];
            dcdecode(temp);
            if(i!=0)
                temp.block[0][0]+=reduct[i-1].block[0][0];
            acdecode(temp);
            reduct.push_back(temp);
        }


        for(int i=0;i<4096;i++)
        {
            deQuantize(reduct[i].block);
            idct2(reduct[i].block, 8);
        }

        FILE *f0 = fopen(outfilename.c_str(), "wb");
        unsigned char enter=0;
        for(int i=0;i<64;i++)
        {
            for(int a=0;a<8;a++)
            {
                for(int j=0;j<64;j++)
                {
                    for(int b=0;b<8;b++)
                    {
                        int put=reduct[i*64+j].block[a][b]+128;
                        enter=(int)put;
                        fwrite(&enter, sizeof(char), 1, f0);
                    }
                }
            }
        }
    }
    return 0;
}
