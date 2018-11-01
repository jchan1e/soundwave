#include <cstdlib>
#include <cmath>
#include <iostream>
//#include <cstring>
#include <fftw3.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <mutex>

using namespace std;

int W = 1920;   // window width
int S = 1080;   // window height
//const int N = W/Q;   // array size



class myRec : public sf::SoundRecorder
{
   struct sramples
   {
      sf::Int16 *samples = NULL;
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

      if (changed)
      {
         sf::Int16* temp = spamples.samples;
         spamples.sampleCount += sampleCount;
         spamples.samples = (sf::Int16*)malloc(spamples.sampleCount*sizeof(sf::Int16));
         int i = 0;
         for (; i < spamples.sampleCount-int(sampleCount); ++i)
            spamples.samples[i] = temp[i];
         for (int j=0; i < spamples.sampleCount; ++i, ++j)
            spamples.samples[i] = samples[j];
         free(temp);
      }

      else
      {
         if (spamples.samples != NULL)
            free(spamples.samples);
         spamples.samples = (sf::Int16*)malloc(sampleCount*sizeof(sf::Int16));
         spamples.sampleCount = (int)sampleCount;
         for (int i=0; i < spamples.sampleCount; ++i)
            spamples.samples[i] = samples[i];
      }

      changed = true;
      m.unlock();

      return true;
   }

   int getSpamples(sf::Int16 array[], int n)
   {
      m.lock();
      for (int i = 0; i < n; ++i)
      {
         if (i < n - spamples.sampleCount)
            array[i] = array[i + spamples.sampleCount];
         else
            array[i] = spamples.samples[spamples.sampleCount - n + i];
      }
      changed = false;
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
void draw_array(int n, double arr[], double bot[]);
void draw_thirds(int n, double arr[]);
int pixel(double wx, double w0, double w1, int p0, int p1);
float lerp(float x1, float x2, float x);
double catmull_rom(double a_y, double b_y, double c_y, double d_y, double t);

//sf::RenderWindow* window;
SDL_Window* window = NULL;
SDL_GLContext context;

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
      for (unsigned int i = 0; i < strang.size(); ++i)
         cout << i << "\t" << strang[i] << endl;
      //cout << rec.getDefaultDevice() << endl << rec.getDevice() << endl;
      int choice = 0;
      if (strang.size() > 1)
         cin >> choice;
      rec.setDevice(rec.getAvailableDevices()[choice]);

      SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER);
      window = SDL_CreateWindow("Soundwave", 0,0, W,S, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
      context = SDL_GL_CreateContext(window);
      SDL_GL_SetSwapInterval(1);


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
      plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_MEASURE);

      int samplesSize = n;
      sf::Int16 samples[samplesSize];
      for (int i = 0; i < samplesSize; ++i)
         samples[i] = 0;

      int i = 0;

      bool quit = false;
      while (!quit)
      {
         SDL_Event event;
         while (SDL_PollEvent(&event))
         {
            if (event.type == SDL_QUIT)
            {
               quit = true;
            }
            else if (event.type == SDL_KEYDOWN && event.key.keysym.scancode == SDL_SCANCODE_Q)
            {
               quit = true;
            }
            else if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
            {
               W = event.window.data1;
               S = event.window.data2;
               f = W;
               free(display);
               display = (double*)malloc(f*sizeof(double));
               glViewport(0,0, W,S);
               glMatrixMode(GL_PROJECTION);
               glLoadIdentity();
               glOrtho(0,W, S,0, -1,1);
               glMatrixMode(GL_MODELVIEW);
               glLoadIdentity();
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

         double f0 = 1.0/(n/88200.0);
         double thirds[33] = {0.0};
         for (i = 1; i < n/2; ++i)
         {
            out[i][0] *= 2.0/n;
            out[i][1] *= 2.0/n;
            out[n-i-1][0] *= 2.0/n;
            out[n-i-1][1] *= 2.0/n;
            processed[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
            processed[i] += out[n-i-1][0]*out[n-i-1][0] + out[n-i-1][1]*out[n-i-1][1];
            //processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? (0.1*((float)i - pow(M_E, 3.0)) + 1) : log(pow((float)i, 1.0/3.0)));
            processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? 1 : log(pow((float)i, 1.0/2.0)));
//            processed[i] = 10.0/log(10.0) * log(processed[i] + 1e-5);
            if (processed[i] < 0.0)
               processed[i] = 0.0;
//            else if (processed[i] > 96.0)
//               processed[i] = 96.0;

            //double fm = n/2.0;
            //if (i == n/2-1)
            //  cout << i << "\t" << fm << endl;
            int ith;
            double ii = f0*i;
            if (17780 < ii && ii <= 22390)
              ith = 32;
            else if (14130 < ii)
              ith = 31;
            else if (11220 < ii)
              ith = 30;
            else if (8913 < ii)
              ith = 29;
            else if (7079 < ii)
              ith = 28;
            else if (5623 < ii)
              ith = 27;
            else if (4467 < ii)
              ith = 26;
            else if (3548 < ii)
              ith = 25;
            else if (2818 < ii)
              ith = 24;
            else if (2239 < ii)
              ith = 23;
            else if (1778 < ii)
              ith = 22;
            else if (1413 < ii)
              ith = 21;
            else if (1122 < ii)
              ith = 20;
            else if (891 < ii)
              ith = 19;
            else if (708 < ii)
              ith = 18;
            else if (562 < ii)
              ith = 17;
            else if (447 < ii)
              ith = 16;
            else if (355 < ii)
              ith = 15;
            else if (282 < ii)
              ith = 14;
            else if (224 < ii)
              ith = 13;
            else if (178 < ii)
              ith = 12;
            else if (141 < ii)
              ith = 11;
            else if (112 < ii)
              ith = 10;
            else if (89.1 < ii)
              ith = 9;
            else if (70.8 < ii)
              ith = 8;
            else if (56.2 < ii)
              ith = 7;
            else if (44.7 < ii)
              ith = 6;
            else if (35.5 < ii)
              ith = 5;
            else if (28.2 < ii)
              ith = 4;
            else if (22.4 < ii)
              ith = 3;
            else if (17.8 < ii)
              ith = 2;
            else if (14.1 < ii)
              ith = 1;
            else // (11.2 < ii)
              ith = 0;

            if (11.2 < ii && ii <= 22390)
              thirds[ith] += processed[i];
         }

         double tfreq[33] = {14.1,17.8,22.4,28.2,35.5,44.7,56.2,70.8,89.1,112,141,178,224,282,355,447,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220,14130,17780};

         double dispth[f];
         int j = 0;
         for (int i=0; i < f; ++i)
         {
            int fi = f0*pow(20000.0/f0, (double)i/f);
            while (j < 32 && fi > tfreq[j])
              ++j;
            dispth[i] = thirds[j];
            //cout << j << endl;
         }

         double A = f0/(rate/n);
         double r = pow(20000.0/(rate/n)/A, 1.0/f);
         //double A = f0;
         //double r = pow(20000.0/f0, 1.0/f);
         // n = A*r^f
         // r = (n/A)^(1/f)
         // 20000 = 20*r^1920
         int bef;
         int pre;
         int nex;
         int aft;
         int i0 = 0;
         int i1 = 1;
         while ((int)(A*pow(r,i0)) == (int)(A*pow(r,i1))) ++i1;
         for (i = 0; i < f; ++i)
         {
            if (i >= i1)
            {
               i0 = i;
               ++i1;
               while ((int)(A*pow(r,i0)) == (int)(A*pow(r,i1))) ++i1;
            }

            pre = A*pow(r, i0);
            bef = max(pre-1, 0);
            nex = min(pre+1, n-1);
            aft = min(nex+1, n-1);
            //cout << pre << endl;

            if (pre < n)
            {
               //display[i] = processed[pre];
               //display[i] = lerp(processed[pre], processed[nex], (float)(i-i0)/(i1-i0));
               display[i] = max(0.0,catmull_rom(processed[bef], processed[pre], processed[nex], processed[aft], (double)(i-i0)/(i1-i0)));
               //cout << processed[bef] << "\t" << processed[pre] << "\t" << processed[nex] << "\t" << processed[aft] << "\t" << (float)(i-i0)/(i1-i0) << "\t" << display[i] << endl;
            }
         }
      //cout << f << endl;
      //cout << n << endl;

         glClear(GL_COLOR_BUFFER_BIT);
         glMatrixMode(GL_MODELVIEW);
         glLoadIdentity();
         draw_array(f, display, dispth);
         draw_wave(n<f/2?n:f/2, arr);
         glFlush();
         SDL_GL_SwapWindow(window);
      }

      rec.stop();
      free(display);
      fftw_destroy_plan(plan);
      fftw_free(in);
      fftw_free(out);
      SDL_Quit();

      return EXIT_SUCCESS;
   }
   //else //(argc == 2) Run on audio file
   //{

// //     else if (argc == 1)
// //        goto MIC;
   //   srand(time(NULL));

   //   //int array[N] = {0};
   //   int i;

   //   sf::ContextSettings settings;
   //   settings.antialiasingLevel = 0;
   //   window = new sf::RenderWindow(sf::VideoMode(W, S), "SoundWave");
   //   window->setVerticalSyncEnabled(true);
   //   //window->setFramerateLimit(60);

   //   sf::SoundBuffer buffer;
   //   if (!buffer.loadFromFile(argv[1]))
   //   {
   //      printf("Could Not Load Audio File\n");
   //      printf("%s\n", argv[1]);
   //      return 1;
   //   }
   //   printf("%s\n", argv[1]);
   //   sf::Sound sound;
   //   sound.setBuffer(buffer);

   //   int rate = buffer.getSampleRate();
   //   printf("SampleRate: %d\n", rate);
   //   const sf::Int16* samples = buffer.getSamples();

   //   int n = rate/15;
   //   int arr[n];
   //   int f = W;
   //   double processed[n];
   //   double* display = (double*)malloc(f*sizeof(double));
   //   int offset = 0;

   //   double factor = rate/1000.0;

   //   double *in;
   //   fftw_complex *out;
   //   fftw_plan plan;

   //   in = (double*) fftw_malloc(sizeof(double) * n);
   //   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
   //   plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_PATIENT);

   //   //for (i = 0; i < n; ++i)
   //   //   arr[i] = (int)samples[i+offset];

   //   sound.play();

   //   sf::Event event;
   //   while (window->isOpen())
   //   {
   //      while (window->pollEvent(event))
   //      {
   //         if (event.type == sf::Event::Closed)
   //         {
   //            window->close();
   //            break;
   //         }
   //         else if (event.type == sf::Event::Resized)
   //         {
   //            W = event.size.width;
   //            S = event.size.height;
   //            f = W;
   //            free(display);
   //            display = (double*)malloc(f*sizeof(double));
   //            window->setView(sf::View(sf::FloatRect(0, 0, event.size.width, event.size.height)));
   //         }
   //      }

   //      if (sound.getStatus() != sf::Music::Playing)
   //      {
   //         window->close();
   //         break;
   //      }

   //      offset = (int) (factor * sound.getPlayingOffset().asMilliseconds());
   //      //printf("%d\n", offset);
   //      for (i = 0; i < n; ++i)
   //      {
   //         if (i >= 0)
   //            arr[i] = avg(buffer.getChannelCount(), samples, i + offset - n);
   //         else
   //            arr[i] = 0;
   //         in[i] = (double)arr[i];
   //      }

   //      fftw_execute(plan);

   //      for (i = 0; i < n; ++i)
   //      {
   //         out[i][0] *= 2.0/n;
   //         out[i][1] *= 2.0/n;
   //         processed[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
   //         //processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? (0.1*((float)i - pow(M_E, 3.0)) + 1) : log(pow((float)i, 1.0/3.0)));
   //         processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? 1 : log(pow((float)i, 1.0/3.0)));
// //           processed[i] = 10.0/log(10.0) * log(processed[i] + 1e-5);
   //         if (processed[i] < 0.0)
   //            processed[i] = 0.0;
// //           else if (processed[i] > 96.0)
// //              processed[i] = 96.0;
   //      }

   //      double A = 40.0/(rate/n);
   //      double r = pow(20000.0/A*n/rate, 1.0/f);
   //      // n = A*r^f
   //      // 20000 = 20*r^1920
   //      int pre;
   //      int nex;
   //      for (i = 0; i < f; ++i)
   //      {
   //         pre = A*pow(r, i+1);
   //         nex = A*pow(r, i+2);
   //         if (nex == pre) ++nex;
   //         display[i] = 0;
   //         //cout << i << ",\t" << pre << endl;
   //         for (int j = pre; j < nex; ++j)
   //         {
   //            if (j < n)
   //               display[i] = max(display[i], processed[j]);
   //         }
   //      }
   //   //cout << f << endl;
   //   //cout << n << endl;

   //      window->clear(sf::Color(0, 0, 0, 12));
   //      draw_wave(n<f/2?n:f/2, arr);
   //      draw_array(f, display);
   //      window->display();
   //   }

   //   free(display);
   //   fftw_destroy_plan(plan);
   //   fftw_free(in);
   //   fftw_free(out);
   //   delete window;

   //   return EXIT_SUCCESS;
   //}
   return EXIT_SUCCESS;
}



//////////// FUNCTION DEFINITIONS /////////////

float lerp(float x0, float x1, float x)
{
  return x0 + (x1-x0)*x;
}

double catmull_rom(double a_y, double b_y, double c_y, double d_y, double t)
{
   return 0.5*((b_y*2.0)+(-a_y+c_y)*t+((a_y*2.0)-(b_y*5.0)+(c_y*4.0)-d_y)*(t*t)+(-a_y+(b_y*3.0)-(c_y*3.0)+d_y)*(t*t*t));
}

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

   glBegin(GL_LINE_STRIP);
   glColor3f(1,1,1);
   //int u = 98304 / S;
   float u = (S/3.0)/(256.0*128.0);
   for(i = 0; i < n; ++i)
   {
      glVertex2f(2*i, S-S/3.0-(arr[i]*u));
   //   line0[1].color = sf::Color::Color (
   //      min(abs(max(abs(((2*arr[i]+512) %1536)-768)-256,0)),255),
   //      min(abs(max(abs(((2*arr[i]) %1536)-768)-256,0)),255),
   //      min(abs(max(abs(((2*arr[i]-512) %1536)-768)-256,0)),255));
   }
   glEnd();
}

void draw_array(int n, double arr[], double bot[])
{
   int i; //y;

   float u = S/300.0;

   glBegin(GL_LINES);
   for(i = 1; i < n; ++i)
   {
      //if(arr[i] > 10.0)
      //   printf("%f ", arr[i]);
      float color_factor = 1920.0*3.0/4.0/W;
      float R = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/24.0, 1.0), 0.20);
      float G = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+1024),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/24.0, 1.0), 0.20);
      float B = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+512),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/24.0, 1.0), 0.20);
      float r = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i),1536)-768)-256,0.0)),256.0);
      float g = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+1024),1536)-768)-256,0.0)),256.0);
      float b = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+512),1536)-768)-256,0.0)),256.0);
   //   line0[1].color = sf::Color::White;
//      line0[0].color = line0[1].color;
      //line0[0] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3-(arr[i]*u)));
      //line0[0].color = line0[1].color;
      //line0[1] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3+(arr[i]*u/2)));
      //line0[1].color = line0[0].color;
      glColor3f(0,0,0);
      glVertex2f(i, S/3.0-(arr[i]*u));
      glColor3f(R,G,B);
      glVertex2f(i, S/3.0);
      glVertex2f(i, S/3.0);
      glColor3f(r,g,b);
      glVertex2f(i, S/3.0+(bot[i]*u/5));
   }
   glEnd();
}

void draw_thirds(int n, double arr[])
{
   int i; //y;

   float u = S/300.0;

   glBegin(GL_LINES);
   for(i = 1; i < n; ++i)
   {
      int j;
      double f0 = 1.0/(n/88200.0);
      double ii = f0*i;
      if (17780 < ii && ii <= 22390)
        j = 32;
      else if (14130 < ii)
        j = 31;
      else if (11220 < ii)
        j = 30;
      else if (8913 < ii)
        j = 29;
      else if (7079 < ii)
        j = 28;
      else if (5623 < ii)
        j = 27;
      else if (4467 < ii)
        j = 26;
      else if (3548 < ii)
        j = 25;
      else if (2818 < ii)
        j = 24;
      else if (2239 < ii)
        j = 23;
      else if (1778 < ii)
        j = 22;
      else if (1413 < ii)
        j = 21;
      else if (1122 < ii)
        j = 20;
      else if (891 < ii)
        j = 19;
      else if (708 < ii)
        j = 18;
      else if (562 < ii)
        j = 17;
      else if (447 < ii)
        j = 16;
      else if (355 < ii)
        j = 15;
      else if (282 < ii)
        j = 14;
      else if (224 < ii)
        j = 13;
      else if (178 < ii)
        j = 12;
      else if (141 < ii)
        j = 11;
      else if (112 < ii)
        j = 10;
      else if (89.1 < ii)
        j = 9;
      else if (70.8 < ii)
        j = 8;
      else if (56.2 < ii)
        j = 7;
      else if (44.7 < ii)
        j = 6;
      else if (35.5 < ii)
        j = 5;
      else if (28.2 < ii)
        j = 4;
      else if (22.4 < ii)
        j = 3;
      else if (17.8 < ii)
        j = 2;
      else if (14.1 < ii)
        j = 1;
      else // (11.2 < ii)
        j = 0;

      //if(arr[i] > 10.0)
      //   printf("%f ", arr[i]);
      float color_factor = 1920.0*3.0/4.0/W;
      float R = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i),1536)-768)-256,0.0)),256.0)*max(min(arr[j]/24.0, 1.0), 0.20);
      float G = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+1024),1536)-768)-256,0.0)),256.0)*max(min(arr[j]/24.0, 1.0), 0.20);
      float B = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+512),1536)-768)-256,0.0)),256.0)*max(min(arr[j]/24.0, 1.0), 0.20);
   //   line0[1].color = sf::Color::White;
//      line0[0].color = line0[1].color;
      //line0[0] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3-(arr[i]*u)));
      //line0[0].color = line0[1].color;
      //line0[1] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3+(arr[i]*u/2)));
      //line0[1].color = line0[0].color;
      //glColor3f(0,0,0);
      //glVertex2f(i, S/3.0-(arr[i]*u));
      glColor3f(R,G,B);
      //glVertex2f(i, S/3.0);
      glVertex2f(i, S/3.0);
      //glColor3f(0,0,0);
      glVertex2f(i, S/3.0+(arr[j]*u/2));
   }
   glEnd();
}

int pixel(double wx, double w0, double w1, int p0, int p1)
{
   double x = (wx-w0)/(w1-w0);
   int px = int (x*(p1-p0)+p0);
   return px;
}
