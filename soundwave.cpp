#include <cstdlib>
#include <cmath>
#include <iostream>
//#include <cstring>
#include <fftw3.h>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <mutex>

using namespace std;

int W = 1920;   // window width
int S = 1080;   // window height
const int Q = 1;   // bar width
//const int N = W/Q;   // array size



class myRec : public sf::SoundRecorder
{
   struct sramples
   {
      sf::Int16 *samples;
      int sampleCount;
   };
   sramples spamples;
   bool changed;
   mutex m;
   
public:

   virtual bool onStart()
   {
      changed = false;
      setProcessingInterval(sf::milliseconds(16));
      return true;
   }

   virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
   {
      m.lock();
      spamples.samples = (sf::Int16*) samples;
      spamples.sampleCount = (int)sampleCount;
      changed = true;
      m.unlock();

      return true;
   }

   int getSpamples(sf::Int16 array[], int n)
   {
      m.lock();
      if (changed)
      {
         for (int i = 0; i < n; ++i)
         {
            if (i < n - spamples.sampleCount)
               array[i] = array[i + spamples.sampleCount];
            else
               array[i] = spamples.samples[spamples.sampleCount - n + i];
         }
         changed = false;
      }
      m.unlock();
      return (spamples.sampleCount);
   }

   bool getChanged()
   {
      return changed;
   }

   void setChanged(bool val)
   {
      changed = val;
   }
};



int avg(int n, const sf::Int16 arr[], int offset);
void draw_wave(int n, int arr[]);
void draw_array(int n, double arr[]);
int pixel(double wx, double w0, double w1, int p0, int p1);

sf::RenderWindow* window;

int main(int argc, char* argv[])
{
   if (argc > 2)
   {
      // Cannot call main as a function.

      //char* argd[1];
      //for (int i = 1; i < argc; ++i)
      //{
      //   argd[0] = argv[i];
      //   main(2, argd);
      //}

      printf("Usage: ./soundwave [audiofile]\n");
      return 1;
   }
   else if (argc < 2) //run on raw mic input
   {
      myRec rec;
      if (!rec.isAvailable())
      {
         printf("System does not support audio capture\n");
         return 2;
      }

      vector<string> strang = rec.getAvailableDevices();
      for (int i = 0; i < strang.size(); ++i)
         cout << strang[i] << endl;
      //cout << rec.getDefaultDevice() << endl << rec.getDevice() << endl;
      int choice = 0;
      if (strang.size() > 1)
         cin >> choice;
      rec.setDevice(rec.getAvailableDevices()[choice]);

      window = new sf::RenderWindow(sf::VideoMode(W, S), "SoundWave");

      sf::ContextSettings settings;
      settings.antialiasingLevel = 0;
      window->setVerticalSyncEnabled(true);

      window->clear(sf::Color(0, 0, 0));
      window->display();
      window->clear(sf::Color(0, 0, 0));
      window->display();

      rec.start(88200);

      int rate = rec.getSampleRate();
      printf("SampleRate: %d\n", rate);
      //const sf::Int16* samples = buffer.getSamples();

      int n = rate/15;
      int f = W;
      int arr[n];
      double processed[n];
      double* display = (double*)malloc(f*sizeof(double));

      //double factor = rate/1000.0;

      double *in;
      fftw_complex *out;
      fftw_plan plan;

      in = (double*) fftw_malloc(sizeof(double) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_PATIENT);
      
      int samplesSize = n;
      sf::Int16 samples[samplesSize];
      for (int i = 0; i < samplesSize; ++i)
         samples[i] = 0;

      int i = 0;

      //rec.start();
      while (window->isOpen())
      {
         sf::Event event;
         while (window->pollEvent(event))
         {
            if (event.type == sf::Event::Closed)
            {
               rec.stop();
               window->close();
               fftw_destroy_plan(plan);
               fftw_free(in);
               fftw_free(out);

               return EXIT_SUCCESS;
            }
            else if (event.type == sf::Event::Resized)
            {
               W = event.size.width;
               S = event.size.height;
               f = W;
               free(display);
               display = (double*)malloc(f*sizeof(double));
               window->setView(sf::View(sf::FloatRect(0, 0, event.size.width, event.size.height)));
            }
         }


         if (rec.getChanged())
         {
            rec.getSpamples(samples, samplesSize);
         }
         //else
         //   offset += n/4;
         //offset = (int) (factor * sound.getPlayingOffset().asMilliseconds());
         //printf("%d\n", offset);
         for (i = 0; i < n; ++i)
         {
            arr[i] = (int)samples[i];
            in[i] = (double)arr[i];
         }
         //cout << "Current: " << current << "\nOffset: " << offset << "\t" << offset + n << endl;

         fftw_execute(plan);

         for (i = 0; i < n; ++i)
         {
            out[i][0] *= 2.0/n;
            out[i][1] *= 2.0/n;
            processed[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
            //processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? (0.1*((float)i - pow(M_E, 3.0)) + 1) : log(pow((float)i, 1.0/3.0)));
            processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? 1 : log(pow((float)i, 1.0/2.0)));
//            processed[i] = 10.0/log(10.0) * log(processed[i] + 1e-5);
            if (processed[i] < 0.0)
               processed[i] = 0.0;
//            else if (processed[i] > 96.0)
//               processed[i] = 96.0;
         }

         double A = 40.0/(rate/n);
         double r = pow(20000.0/A*n/rate, 1.0/f);
         // n = A*r^f
         // r = (n/A)^(1/f)
         // 20000 = 20*r^1920
         int pre;
         int nex;
         for (i = 0; i < f; ++i)
         {
            pre = A*pow(r, i+1);
            nex = A*pow(r, i+2);
            if (nex == pre) ++nex;
            display[i] = 0;
            //cout << i << ",\t" << pre << endl;
            for (int j = pre; j < nex; ++j)
            {
               if (j < n)
                  display[i] = max(display[i], processed[j]);
            }
         }
      //cout << f << endl;
      //cout << n << endl;

         window->clear(sf::Color(0, 0, 0));
         
         draw_wave((n-1)<f/2?(n-1):f/2, arr);
         draw_array(f, display);
         window->display();
         //draw_array(f, display);
         //draw_wave(n, arr);
         //window->display();
      }

      free(display);
      fftw_destroy_plan(plan);
      fftw_free(in);
      fftw_free(out);
      delete window;

      return EXIT_SUCCESS;
   }
   else //(argc == 2) Run on audio file
   {

//      else if (argc == 1)
//         goto MIC;
      srand(time(NULL));

      //int array[N] = {0};
      int i;

      sf::ContextSettings settings;
      settings.antialiasingLevel = 0;
      window = new sf::RenderWindow(sf::VideoMode(W, S), "SoundWave");
      window->setFramerateLimit(60);

      sf::SoundBuffer buffer;
      if (!buffer.loadFromFile(argv[1]))
      {
         printf("Could Not Load Audio File\n");
         printf("%s\n", argv[1]);
         return 1;
      }
      printf("%s\n", argv[1]);
      sf::Sound sound;
      sound.setBuffer(buffer);

      int rate = buffer.getSampleRate();
      printf("SampleRate: %d\n", rate);
      const sf::Int16* samples = buffer.getSamples();

      int n = rate/15;
      int arr[n];
      int f = W;
      double processed[n];
      double* display = (double*)malloc(f*sizeof(double));
      int offset = 0;

      double factor = rate/1000.0;

      double *in;
      fftw_complex *out;
      fftw_plan plan;

      in = (double*) fftw_malloc(sizeof(double) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_PATIENT);

      //for (i = 0; i < n; ++i)
      //   arr[i] = (int)samples[i+offset];

      sound.play();

      sf::Event event;
      while (window->isOpen())
      {
         while (window->pollEvent(event))
         {
            if (event.type == sf::Event::Closed)
            {
               window->close();
               break;
            }
            else if (event.type == sf::Event::Resized)
            {
               W = event.size.width;
               S = event.size.height;
               f = W;
               free(display);
               display = (double*)malloc(f*sizeof(double));
               window->setView(sf::View(sf::FloatRect(0, 0, event.size.width, event.size.height)));
            }
         }

         if (sound.getStatus() != sf::Music::Playing)
         {
            window->close();
            break;
         }

         offset = (int) (factor * sound.getPlayingOffset().asMilliseconds());
         //printf("%d\n", offset);
         for (i = 0; i < n; ++i)
         {
            if (i >= 0)
               arr[i] = avg(buffer.getChannelCount(), samples, i + offset - n);
            else
               arr[i] = 0;
            in[i] = (double)arr[i];
         }

         fftw_execute(plan);

         for (i = 0; i < n; ++i)
         {
            out[i][0] *= 2.0/n;
            out[i][1] *= 2.0/n;
            processed[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
            //processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? (0.1*((float)i - pow(M_E, 3.0)) + 1) : log(pow((float)i, 1.0/3.0)));
            processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? 1 : log(pow((float)i, 1.0/3.0)));
//            processed[i] = 10.0/log(10.0) * log(processed[i] + 1e-5);
            if (processed[i] < 0.0)
               processed[i] = 0.0;
//            else if (processed[i] > 96.0)
//               processed[i] = 96.0;
         }

         double A = 40.0/(rate/n);
         double r = pow(20000.0/A*n/rate, 1.0/f);
         // n = A*r^f
         // 20000 = 20*r^1920
         int pre;
         int nex;
         for (i = 0; i < f; ++i)
         {
            pre = A*pow(r, i+1);
            nex = A*pow(r, i+2);
            if (nex == pre) ++nex;
            display[i] = 0;
            //cout << i << ",\t" << pre << endl;
            for (int j = pre; j < nex; ++j)
            {
               if (j < n)
                  display[i] = max(display[i], processed[j]);
            }
         }
      //cout << f << endl;
      //cout << n << endl;

         window->clear(sf::Color(0, 0, 0, 12));
         draw_wave((n-1)<f/2?(n-1):f/2, arr);
         draw_array(f, display);
         window->display();
      }

      free(display);
      fftw_destroy_plan(plan);
      fftw_free(in);
      fftw_free(out);
      delete window;

      return EXIT_SUCCESS;
   }
}



//////////// FUNCTION DEFINITIONS /////////////


int avg(int n, const sf::Int16 arr[], int offset)
{
   int a = 0;
   for(int i = 0; i < n; ++i)
      a += arr[i + n*offset];
   return a/n;
}

void draw_wave(int n, int arr[])
{
   int i; //y;

   sf::Vertex* line0 = (sf::Vertex*)malloc(2*n*sizeof(sf::Vertex));
   //{  
   //   sf::Vertex(sf::Vector2f(0, 0)),
   //   sf::Vertex(sf::Vector2f(0, 0)),
   //};

   //int u = 98304 / S;
   float u = (S/3.0)/(256.0*128.0);
   for(i = 0; i < n; ++i)
   {
      //line0[2*i] = sf::Vertex(sf::Vector2f(2*(i-1), S-S/3.0));
      //line0[2*i+1] = sf::Vertex(sf::Vector2f(2*i, S-S/3.0));
      line0[2*i] = sf::Vertex(sf::Vector2f(2*i, S-S/3.0-(arr[i]*u)));
      line0[2*i+1] = sf::Vertex(sf::Vector2f(2*(i+1), S-S/3.0-(arr[i+1]*u)));
   //   line0[1].color = sf::Color::Color (
   //      min(abs(max(abs(((2*arr[i]+512) %1536)-768)-256,0)),255),
   //      min(abs(max(abs(((2*arr[i]) %1536)-768)-256,0)),255),
   //      min(abs(max(abs(((2*arr[i]-512) %1536)-768)-256,0)),255));
      line0[2*i].color = sf::Color::White;
      line0[2*i+1].color = sf::Color::White;
//      line0[0].color = sf::Color::Black;
   }
   window->draw(line0, 2*n, sf::Lines);//Strip);
   free(line0);
}

void draw_array(int n, double arr[])
{
   int i; //y;

   sf::Vertex line0[4];
   //{  
   //   sf::Vertex(sf::Vector2f(0, 0)),
   //   sf::Vertex(sf::Vector2f(0, 0)),
   //};

   float u = S/300.0;

   for(i = 1; i < n; ++i)
   {
      //if(arr[i] > 10.0)
      //   printf("%f ", arr[i]);
      for(int j = 0; j < Q; ++j)
      {
         line0[0] = sf::Vertex(sf::Vector2f(Q*i+j, S/3.0-(arr[i]*u)));
         line0[1] = sf::Vertex(sf::Vector2f(Q*i+j, S/3.0));
         line0[2] = sf::Vertex(sf::Vector2f(Q*i+j, S/3.0));
         line0[3] = sf::Vertex(sf::Vector2f(Q*i+j, S/3.0+(arr[i]*u/2)));
         line0[1].color = sf::Color (
            min(abs(max(abs(((2*i/3) %1536)-768)-256,0)),255)*max(min(arr[i]/24.0, 1.0), 0.20),
            min(abs(max(abs(((2*i/3+1024) %1536)-768)-256,0)),255)*max(min(arr[i]/24.0, 1.0), 0.20),
            min(abs(max(abs(((2*i/3+512) %1536)-768)-256,0)),255)*max(min(arr[i]/24.0, 1.0), 0.20));
   //      line0[1].color = sf::Color::White;
//         line0[0].color = line0[1].color;
         line0[2].color = line0[1].color;
         line0[0].color = sf::Color(0,0,0,0);
         line0[3].color = sf::Color(0,0,0,0);
         //window->draw(line0, 2, sf::Lines);//Strip);
         //line0[0] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3-(arr[i]*u)));
         //line0[0].color = line0[1].color;
         //line0[1] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3+(arr[i]*u/2)));
         //line0[1].color = line0[0].color;
         window->draw(line0, 4, sf::Lines);
      }
      //window->draw(line0, 2*n, sf::Lines);//Strip);
   }
}

int pixel(double wx, double w0, double w1, int p0, int p1)
{
   double x = (wx-w0)/(w1-w0);
   int px = int (x*(p1-p0)+p0);
   return px;
}
