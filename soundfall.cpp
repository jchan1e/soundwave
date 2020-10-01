#include <cstdlib>
#include <cmath>
#include <iostream>
#include <deque>
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
unsigned int frame;


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
void draw_array(int n, double arr[]);
void draw_map();
void increment_map(int n, double* row);
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
      window = SDL_CreateWindow("Soundfall", 0,0, W,S, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
      context = SDL_GL_CreateContext(window);
      if (SDL_GL_SetSwapInterval(1) < 0)
        cout << "SDL could not set Vsync: " << SDL_GetError() << endl;

      glGenTextures(1,&frame);
      glBindTexture(GL_TEXTURE_2D,frame);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);

      // Initialize map texture to black
      glClear(GL_COLOR_BUFFER_BIT);
      glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,0,0,W,S,0);

      rec.start(88200);

      int rate = rec.getSampleRate();
      printf("SampleRate: %d\n", rate);
      //const sf::Int16* samples = buffer.getSamples();

      int n = rate/15;
      int f = W;
      int arr[n];
      double processed[n];
      double* display;
      //unsigned int v_size = S - S/3;
      //deque<float*> vec;

      display = (double*)malloc(f*sizeof(double));

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
               //v_size = S-S/3;
               free(display);
               display = (double*)malloc(f*sizeof(double));
               //for (float* v : vec) {
               //  free(v);
               //  vec.pop_front();
               //}
               glViewport(0,0, W,S);
               glMatrixMode(GL_PROJECTION);
               glLoadIdentity();
               glOrtho(0,W, S,0, -1,1);
               glMatrixMode(GL_MODELVIEW);
               glLoadIdentity();

               glDeleteTextures(1,&frame);
               glGenTextures(1,&frame);
               glBindTexture(GL_TEXTURE_2D,frame);
               glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
               glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
               glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
               glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);

               glClear(GL_COLOR_BUFFER_BIT);
               glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,0,0,W,S,0);
            }
         }

         if (rec.getChanged())
         {
            rec.getSpamples(samples, samplesSize);
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
               processed[i] = powf(out[i][0]*out[i][0] + out[i][1]*out[i][1], 0.5);
               //processed[i] = powf(processed[i], 1.0/2.0)/60 * (i<pow(M_E, 3.0) ? (0.1*((float)i - pow(M_E, 3.0)) + 1) : log(pow((float)i, 1.0/3.0)));
               processed[i] = processed[i] / 240 * sqrt((float)i);
//               processed[i] = 10.0/log(10.0) * log(processed[i] + 1e-5);
               if (processed[i] < 0.0)
                  processed[i] = 0.0;
//               else if (processed[i] > 96.0)
//                  processed[i] = 96.0;
            }

            double A = 40.0/(rate/n);
            double r = pow(20000.0/A*n/rate, 1.0/f);
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
               bef = max(pre-1, 0);;
               nex = pre+1;
               aft = nex+1;

               //display[i] = processed[pre];
               //display[i] = lerp(processed[pre], processed[nex], (float)(i-i0)/(i1-i0));
               display[i] = max(0.0,catmull_rom(processed[bef], processed[pre], processed[nex], processed[aft], (double)(i-i0)/(i1-i0)));
               //cout << processed[bef] << "\t" << processed[pre] << "\t" << processed[nex] << "\t" << processed[aft] << "\t" << (float)(i-i0)/(i1-i0) << "\t" << display[i] << endl;
            }
            increment_map(f, display);
            //cout << vec.size() << endl;
         }
      //cout << f << endl;
      //cout << n << endl;

         //vec.push_front(display);
         //if (vec.size() > v_size) {
         //   free(vec[vec.size()-1]);
         //   vec.pop_back();
         //}

         glClear(GL_COLOR_BUFFER_BIT);
         glMatrixMode(GL_MODELVIEW);
         glLoadIdentity();
         draw_map();
         draw_array(f, display);
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

void draw_array(int n, double arr[])
{
   int i; //y;

   float u = S/300.0;

   glBegin(GL_LINES);
   for(i = 1; i < n; ++i)
   {
      //if(arr[i] > 10.0)
      //   printf("%f ", arr[i]);
      float color_factor = 1920.0*3.0/4.0/W;
      float R = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i     ),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/48.0, 1.0), 0.15);
      float G = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+1024),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/48.0, 1.0), 0.15);
      float B = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+512 ),1536)-768)-256,0.0)),256.0)*max(min(arr[i]/48.0, 1.0), 0.15);
      if (arr[i] > 48.0) {
         float R1 = R + (1.0-R) * (arr[i]-48.0)*0.5/48.0;
         float G1 = G + (1.0-G) * (arr[i]-48.0)*0.5/48.0;
         float B1 = B + (1.0-B) * (arr[i]-48.0)*0.5/48.0;

         glColor3f(0,0,0);
         glVertex2f(i, S/3.0-(arr[i]*u));
         glColor3f(R,G,B);
         glVertex2f(i, S/3.0-((arr[i]-48.0)*0.333*u));
         glVertex2f(i, S/3.0-((arr[i]-48.0)*0.333*u));
         glColor3f(R1,G1,B1);
         glVertex2f(i, S/3.0);
      }
      else {
   //      line0[1].color = sf::Color::White;
//         line0[0].color = line0[1].color;
         //line0[0] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3-(arr[i]*u)));
         //line0[0].color = line0[1].color;
         //line0[1] = sf::Vertex(sf::Vector2f(2*Q*i+Q+j, S/3+(arr[i]*u/2)));
         //line0[1].color = line0[0].color;
         glColor3f(0,0,0);
         glVertex2f(i, S/3.0-(arr[i]*u));
         glColor3f(R,G,B);
         glVertex2f(i, S/3.0);
         //glVertex2f(i, S/3.0);
         //glColor3f(0,0,0);
         //glVertex2f(i, S/3.0+(arr[i]*u/12));
      }
   }
   glEnd();
}

void increment_map(int n, double* row) {
   float color_factor = 1920.0*3.0/4.0/W;
   float* arr = (float*)malloc(n*3*sizeof(float));
   for (int i=0; i < n; ++i) {
      float R = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i     ),1536)-768)-256,0.0)),256.0)*max(min(row[i]/48.0, 1.0), 0.0);
      float G = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+1024),1536)-768)-256,0.0)),256.0)*max(min(row[i]/48.0, 1.0), 0.0);
      float B = 1.0/256.0*min(abs(max(abs(fmod((color_factor*i+512 ),1536)-768)-256,0.0)),256.0)*max(min(row[i]/48.0, 1.0), 0.0);
      if (row[i] > 48.0) {
        R += (1.0-R) * (row[i]-48.0)*0.5/48.0;
        G += (1.0-G) * (row[i]-48.0)*0.5/48.0;
        B += (1.0-B) * (row[i]-48.0)*0.5/48.0;
      }
      arr[i*3+0] = R;
      arr[i*3+1] = G;
      arr[i*3+2] = B;
   }
   //redraw map with offset
   glClear(GL_COLOR_BUFFER_BIT);
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D,frame);

   glBegin(GL_QUADS);
   glTexCoord2f(0.0, 2.0/3.0);
   glVertex2f(0.0, S/3+1);
   glTexCoord2f(0.0, 1.0/S);
   glVertex2f(0.0, S);
   glTexCoord2f(1.0, 1.0/S);
   glVertex2f(W, S);
   glTexCoord2f(1.0, 2.0/3.0);
   glVertex2f(W, S/3+1);
   glEnd();
   glDisable(GL_TEXTURE_2D);

   //draw array into the gap
   glBegin(GL_POINTS);
   for (int i=0; i < n; ++i) {
     float R = arr[i*3+0];
     float G = arr[i*3+1];
     float B = arr[i*3+2];
     glColor3f(R,G,B);
     glVertex2f(i, S/3+1);
   }
   glEnd();

   //recapture texture
   glFlush();
   glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,0,0,W,S,0);

   free(arr);
}

void draw_map() {
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D,frame);

   glBegin(GL_QUADS);
   glColor3f(1.0,1.0,1.0);
   glTexCoord2f(0.0, 2.0/3.0);
   glVertex2f(0.0, S/3);
   glTexCoord2f(0.0, 0.0);
   glVertex2f(0.0, S);
   glTexCoord2f(1.0, 0.0);
   glVertex2f(W, S);
   glTexCoord2f(1.0, 2.0/3.0);
   glVertex2f(W, S/3);
   glEnd();

   glDisable(GL_TEXTURE_2D);

   //int v_lim = min((int)vec.size(), S-S/3);
   //glBegin(GL_POINTS);
   //for (int v=0; v < v_lim; ++v) {
   //   for (int i=0; i < n; ++i) {
   //      float R = vec[v][i*3+0];
   //      float G = vec[v][i*3+1];
   //      float B = vec[v][i*3+2];
   //      glColor3f(R,G,B);
   //      glVertex2f(i, v+S/3);
   //   }
   //}
   //glEnd();
}

int pixel(double wx, double w0, double w1, int p0, int p1)
{
   double x = (wx-w0)/(w1-w0);
   int px = int (x*(p1-p0)+p0);
   return px;
}
