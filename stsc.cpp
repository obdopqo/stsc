#include <cstdio>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <queue>

// ÿ���������(��ȷ������֮һ�룬̫�����ʶ����������ǿ��ܻᶪ������)
#define SPR_CNT 20
// ����������������(����Խ����Խ�࣬���99)
#define HFILTER 40
// ��С����(������)
#define MIN_RES 2
// �ﵽ�����仯ʱ�����趨����(�ڽ�ǿ/�����ĸ������ã������Ұ���ֱ�ӵ���100��Ҳ���ǲ������趨)
#define MIN_JMP 100

using namespace std;

const double pi = 3.14159265358979323846;

//fft��������https://blog.csdn.net/enjoy_pascal/article/details/81478582
typedef complex<double> cp;
void fft(cp *a,long n,long inv) //n������2�������� 
{
    long bit=0;
    while ((1<<bit)<n) bit++;
    unsigned long rev[1<<bit] = {};
    for (unsigned long i=0;i<n;i++)
    {
        rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
        if (i<rev[i])swap(a[i],a[rev[i]]);//��������if�ύ�����Σ�����û������
    }
    for (long mid=1;mid<n;mid*=2)//mid��׼���ϲ����еĳ��ȵĶ���֮һ
    {
    	cp temp(cos(pi/mid),inv*sin(pi/mid));//��λ����pi��ϵ��2�Ѿ�Լ����
        for (long i=0;i<n;i+=mid*2)//mid*2��׼���ϲ����еĳ��ȣ�i�Ǻϲ�������һλ
		{
            cp omega(1,0);
            for (long j=0;j<mid;j++,omega*=temp)//ֻɨ��벿�֣��õ��Ұ벿�ֵĴ�
            {
                cp x=a[i+j],y=omega*a[i+j+mid];
                a[i+j]=x+y,a[i+j+mid]=x-y;//������Ǻ����任ʲô��
            }
        }
    }
}

class outnote {
	public:
	outnote(double a,int b,int c,double d):start(a),pitch(b),volume(c),time(d){};
	double start;
	int pitch;
	int volume;
	double time;
};

class mycompare {
	public:
	bool operator()(const outnote &a,const outnote &b){
		return a.start > b.start;
	}
};

int main(int argc,char* argv[]){
	struct riff {
		long ID;
		long Size;
		long Type;
	} riff;
	struct format {
		long ID;
		long Size;
		short AF;//AudioFormat (1)
		short NCs;//NumChannels (1)
		long SR;//SampleRate (11025)
		long BR;//ByteRate (?)
		short BA;//BlockAlign (?)
		short BPS;//BitsPerSample (8)
	}format;
	struct data {
		long ID;
		long Size;
	} data;
	unsigned char* wavraw;
	int wavsize;
	int status[131],status2[131];
	double statust[131];
	FILE* fp=fopen(argc==1?"example.wav":argv[1],"rb");
	if(fp==NULL){
		perror("While opening wave file");
		return 1;
	}

	fread(&riff,sizeof(struct riff),1,fp);
	printf("%08X,%08d,%08X\n",riff.ID,riff.Size,riff.Type);
	fread(&format,sizeof(struct format),1,fp);
	printf("%08X,%08d,NCs=%d,SR=%d,BPS=%d\n",format.ID,format.Size,format.NCs,format.SR,format.BPS);
	fread(&data,sizeof(struct data),1,fp);
	printf("%08X,%08d\n",data.ID,data.Size);

	data.Size=riff.Size+4-sizeof(struct riff)-sizeof(struct format)-sizeof(struct data);
	wavsize = data.Size/format.NCs/(format.BPS/8);
	wavraw = new unsigned char[wavsize];
	//fread(wavraw,sizeof(char),data.Size,fp);
	
	for(int i=0;i<wavsize;i++){
		int v=0;
		for(int j=0;j<format.NCs;j++){
			for(int k=1;k<format.BPS/8;k++){
				fgetc(fp);
			}
			v+=fgetc(fp);
		}
		v/=format.NCs;
		wavraw[i]=v;
	}
	fclose(fp);

	for(int i=0;i<131;i++){
		status[i]=0;
		statust[i]=0;
	}

	int units=wavsize * SPR_CNT / format.SR;
	/* --- when it comes to a world of Scratch --- */
	printf("units=%d\n",units);
	//calculate time
	
	//each time step

	int sample=2;
	while(sample<format.SR/SPR_CNT){
		sample*=2;
	}
	//sample/=2;
	printf("sample=%d\n",sample);

	cp a[sample]={};
	
	FILE* wo=fopen("stsc.txt","w");
		
	priority_queue<outnote,vector<outnote>,mycompare> notes;
	
	for(int tp=0;tp<units-1;tp++){
		int pos=tp*format.SR/SPR_CNT;
		printf("tp=%08d pos=%08d size=%08d proc=%05.2lf%%\r",tp,pos,wavsize,100.0*tp/units-1);
		for(int i=0;i<131;i++){
			status2[i]=0;
		}
		for(int i=0;i<sample;i++){
			a[i].real(wavraw[i+pos]);
			a[i].imag(0);
		}
		fft(a,sample,1);
		for(int i=1;i<sample/2;i++){
			double feq=(double)(i*format.SR)/sample;
			int note=floor(60+log(feq/261.63)/log(2)*12+0.5);
			int res=floor(sqrt(a[i].real()*a[i].real()+a[i].imag()*a[i].imag())/sample*100/255);
			//printf(" feq=%6.2f note=%d res=%d\n",feq,note,res);
			if(note>=0&&note<=130){
				if(res>status2[note]){
					status2[note]=res;
				}
			}
		}
		int high=0;
		for(int note=0;note<131;note++){
			if(status2[note]>high){
				high=status2[note];
			}
		}
		for(int note=0;note<131;note++){
			if(status2[note]<max(high*HFILTER/100,MIN_RES)){
				status2[note]=0;
			}
			if(abs(status2[note]-status[note])>MIN_JMP||status2[note]==0){
				if(status[note]!=0){
					notes.push(outnote(
						statust[note],
						status[note],
						note,
						(double)tp/SPR_CNT-statust[note]
					));
				}
				status[note]=status2[note];
				statust[note]=(double)tp/SPR_CNT;
			}else{
				if(status[note]<status2[note]){
					status[note]=status2[note];
				}
			}
		}
	}

	for(int note=0;note<131;note++){
		if(status[note]!=0){
			notes.push(outnote(
				statust[note],
				status[note],
				note,
				(double)(units-1)/SPR_CNT-statust[note]
			));
		}
	}
	
	while(!notes.empty()){
		const outnote &n=notes.top();
		fprintf(wo,"%lf\n",n.start);
		fprintf(wo,"%d\n",n.pitch);
		fprintf(wo,"%d\n",n.volume);
		fprintf(wo,"%lf\n",n.time);
		notes.pop();
	}

	fclose(wo);
	
	return 0;
}
