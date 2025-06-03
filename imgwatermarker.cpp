// SDL2 + fftw example for DIP2020 course
// 2020 (c) A.Łukaszewski, use as you like
// Compile:
// g++ -o sdl2-fft sdl2-fft.cpp -Wall -lSDL2_image -lSDL2 -lfftw3

#include <iostream>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include <math.h>
#include <fftw3.h>
// Auxilary FFT functions  ----------------------------------------
typedef struct
{
    double re, im;
} zomplex;

fftw_plan *plan;

void nosgi_fft(zomplex *array, int width, int height, int inv)
{
    fftw_complex *carray = (fftw_complex *)malloc(width * height * sizeof(fftw_complex));
    int ind;

    for (ind = 0; ind < width * height; ind++)
    {
        carray[ind][0] = array[ind].re;
        carray[ind][1] = array[ind].im;
    }

    plan[-1] = fftw_plan_dft_2d(height, width, carray, carray, FFTW_BACKWARD, //(fftw_direction)-1,
                                FFTW_ESTIMATE);                               //| FFTW_IN_PLACE);
    plan[1] = fftw_plan_dft_2d(height, width, carray, carray, FFTW_FORWARD,   //(fftw_direction)1,
                               FFTW_ESTIMATE);                                // | FFTW_IN_PLACE);
    fftw_execute(plan[inv]);

    for (ind = 0; ind < width * height; ind++)
    {
        array[ind].re = carray[ind][0];
        array[ind].im = carray[ind][1];
    }

    free(carray);
}

void compute_fft(zomplex *array, int width, int height)
{
    // void initialise_fft (int width, int height) {
    plan = ((fftw_plan *)malloc(3 * sizeof(fftw_plan))) + 1;

    nosgi_fft(array, width, height, -1);
}

void compute_inverse_fft(zomplex *array, int width, int height)
{
    int ind;
    int value = width * height;

    nosgi_fft(array, width, height, 1);
    for (ind = 0; ind < value; ind++)
        array[ind].re /= value;
}

void gaussian_filter(zomplex *filter, int width, int height, double scale)
{
    static double k = 1. / (2. * 1.4142136);
    int x, y, i;
    double x1, y1, s;
    double a = 1. / (k * scale);
    double c = 1. / 4.;

    for (i = 0, y = 0; y < height; y++)
    {
        y1 = (y >= height / 2) ? y - height : y;
        s = erf(a * (y1 - .5)) - erf(a * (y1 + .5));
        for (x = 0; x < width; x++, i++)
        {
            x1 = (x >= width / 2) ? x - width : x;
            filter[i].re = s * (erf(a * (x1 - .5)) - erf(a * (x1 + .5))) * c;
            filter[i].im = 0.;
        }
    }
}
//------------------------------------------------------------------------

int MIN(int x, int y)
{
    return (x < y ? x : y);
}

unsigned char *ImgPtr(int x, int y, SDL_Surface *s)
{ // For 1byte/channel 3-4
    return (unsigned char *)(s->pixels) + s->format->BytesPerPixel * (x + y * s->w);
}

void DoFFT(SDL_Surface *img, SDL_Surface *fft)
{
    // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    int width = img->w, height = img->h;
    // int      len=width*height;
    zomplex *my_fft; //, *kernel;
    int x, y;

    my_fft = (zomplex *)calloc((int)width * height, sizeof(zomplex));
    // kernel =  (zomplex*) calloc ((int)width*height, sizeof (zomplex));
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            my_fft[x + width * y].re =
                0.3 * ImgPtr(x, y, img)[0] + 0.6 * ImgPtr(x, y, img)[1] + 0.1 * ImgPtr(x, y, img)[2];
            my_fft[x + width * y].im = 0;
        }

    // initialise_fft(width,height);
    compute_fft(my_fft, width, height);

    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            double r = sqrt(my_fft[x + width * y].re * my_fft[x + width * y].re +
                            my_fft[x + width * y].im * my_fft[x + width * y].im); //,
            // phi=atan2(my_fft[x+width*y].re, my_fft[x+width*y].im);

            ImgPtr(x, y, fft)[0] = ImgPtr(x, y, fft)[1] = ImgPtr(x, y, fft)[2] =
                MIN(255, (int)(1024. * log(1 + r / (width * height))));
        }
    // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
}
//------------------------------------------------------------------------

//=====================
// Helper stuff
//=====================

void setColor(SDL_Surface *surface, int x, int y, unsigned char r, unsigned char g, unsigned char b){
    ImgPtr(x, y, surface)[0] = r;
    ImgPtr(x, y, surface)[1] = g;
    ImgPtr(x, y, surface)[2] = b;
}





//=====================
// Perform operation
//=====================

void drawFFT(zomplex *fft_data, SDL_Surface *output, int width, int height, datadisplayed_t mode)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int index = x + width * y;
            double re = fft_data[index].re;
            double im = fft_data[index].im;
            unsigned char pixel_value = 0;
            double value;

            switch (mode)
            {
            case datadisplayed_t::DATADISPLAYED_REAL:
                value = re;
                pixel_value = MIN(255, (int)(128 + 20 * log(1 + fabs(value))));
                
                setColor(output, x, y
                    ,value < 0 ? 0 : pixel_value
                    ,value >= 0 ? 0 : pixel_value
                    ,0
                );
                break;
            case datadisplayed_t::DATADISPLAYED_IMAGINARY:
                value = im;
                pixel_value = MIN(255, (int)(128 + 20 * log(1 + fabs(value))));
                
                setColor(output, x, y
                    ,value < 0 ? 0 : pixel_value
                    ,pixel_value
                    ,value >= 0 ? 0 : pixel_value
                );
                break;
            case datadisplayed_t::DATADISPLAYED_MAGNITUDE:
                value = sqrt(re * re + im * im);
                pixel_value = MIN(255, (int)(10 * log(1 + value)));
                
                setColor(output, x, y
                    ,pixel_value
                    ,pixel_value
                    ,pixel_value
                );
                break;
            case datadisplayed_t::DATADISPLAYED_PHASE: // courtesy of some kind stackoverflow folks
            {
                double phase = atan2(im, re);
                // Map phase [-π, π] to [0, 255]
                pixel_value = (unsigned char)(((phase + M_PI) * 128.0) / M_PI);


                // Create a hue-based representation of phase
                double hue = (phase + M_PI) / (2 * M_PI); // Normalize to [0, 1]

                // Simple HSV to RGB conversion (with S=1, V=1)
                double r = 0, g = 0, b = 0;
                int h_i = (int)(6 * hue);
                double f = 6 * hue - h_i;
                switch (h_i % 6)
                {
                case 0:
                    r = 1;
                    g = f;
                    b = 0;
                    break;
                case 1:
                    r = 1 - f;
                    g = 1;
                    b = 0;
                    break;
                case 2:
                    r = 0;
                    g = 1;
                    b = f;
                    break;
                case 3:
                    r = 0;
                    g = 1 - f;
                    b = 1;
                    break;
                case 4:
                    r = f;
                    g = 0;
                    b = 1;
                    break;
                case 5:
                    r = 1;
                    g = 0;
                    b = 1 - f;
                    break;
                }

                setColor(output, x, y
                    ,255 * r
                    ,255 * g
                    ,255 * b
                );                
            }
                break;
            default:
                setColor(output, x, y, 0, 0, 0);
            }
        }
    }
}

void perform_operation
    (SDL_Surface *img
    ,SDL_Surface *fft
    ,SDL_Surface *postOperationFFT
    ,SDL_Surface *postOperationImg
    ,operationtype_t operationType
    ,operationdata_t operationData
    ,datadisplayed_t dataDisplayed){

    int width = img->w, height = img->h;
    // int      len=width*height;
    zomplex *my_fft, *outputFFT, *copy_fft; //, *kernel;
    int x, y;

    my_fft = (zomplex *)calloc((int)width * height, sizeof(zomplex));
    copy_fft = (zomplex *)calloc((int)width * height, sizeof(zomplex));
    outputFFT = (zomplex *)calloc((int)width * height, sizeof(zomplex));



    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            my_fft[x + width * y].re =
                0.3 * ImgPtr(x, y, img)[0] + 0.6 * ImgPtr(x, y, img)[1] + 0.1 * ImgPtr(x, y, img)[2];

            my_fft[x + width * y].im = 0;

            copy_fft[x + width * y].re = my_fft[x + width * y].re;
            copy_fft[x + width * y].im = my_fft[x + width * y].im;

            if(operationType == operationtype_t::ALTERNATING_SIGNS){
                my_fft[x + width * y].re *= 0.5;

                if((x + y) % 2)
                    my_fft[x + width * y].re = 0.5 - my_fft[x + width * y].re;
                else
                    my_fft[x + width * y].re = 0.5 + my_fft[x + width * y].re;
            }
        }

    compute_fft(my_fft, width, height);
    compute_fft(copy_fft, width, height);

    drawFFT(copy_fft, fft, width, height, dataDisplayed);

    // do the operation

    printInfo(operationType, dataDisplayed);
    

    operationData.m_inputfft = my_fft;
    operationData.m_outputfft = outputFFT;
    operationData.m_height = height;
    operationData.m_width = width;

    availableOperations[static_cast<int>(operationType)](operationData);


    drawFFT(outputFFT, postOperationFFT, width, height, dataDisplayed);

    compute_inverse_fft(outputFFT, width, height);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int index = x + width * y;
            double value = outputFFT[index].re;

            // Clamp values to [0, 255]
            unsigned char pixel_value = (unsigned char)std::max(0, std::min(255, (int)value));

            setColor(postOperationImg, x, y, pixel_value, pixel_value, pixel_value);
        }
    }
}



//=====================
// Main program
//=====================

void handle_input(SDL_Event &event, operationtype_t &operationType, datadisplayed_t &dataDisplayed, bool &quit){
    switch (event.key.keysym.sym){
        case SDLK_ESCAPE:
            quit = true;
            break;
        case SDLK_RIGHT:
            {
                int ind = static_cast<int>(dataDisplayed);
                ind = (ind + 1) % static_cast<int>(datadisplayed_t::DATADISPLAYED_COUNT);
                dataDisplayed = static_cast<datadisplayed_t>(ind);
            }
            break;
        case SDLK_LEFT:
            {
                int ind = static_cast<int>(dataDisplayed);
                ind = (ind - 1);
                if(ind < 0) ind += static_cast<int>(datadisplayed_t::DATADISPLAYED_COUNT);
                dataDisplayed = static_cast<datadisplayed_t>(ind);
            }
            break;
        case SDLK_UP:
            {
                int ind = static_cast<int>(operationType);
                ind = (ind + 1) % static_cast<int>(operationtype_t::OPERATIONTYPE_COUNT);
                operationType = static_cast<operationtype_t>(ind);
            }
            break;
        case SDLK_DOWN:
            {
                int ind = static_cast<int>(operationType);
                ind = (ind - 1);
                if(ind < 0) ind += static_cast<int>(operationtype_t::OPERATIONTYPE_COUNT);
                operationType = static_cast<operationtype_t>(ind);
            }
            break;
        case SDLK_0:
            operationType = operationtype_t::DEFAULT;
            break;
        case SDLK_1:
            operationType = operationtype_t::ZERO_IMAGINARY;
            break;
        // ...
    }
}



int main(int argc, char **argv)
{
    bool quit = false;
    SDL_Event event;
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Surface *image = IMG_Load(argc < 2 ? "lena.jpg" : argv[1]);
    SDL_Surface *imfft = IMG_Load(argc < 2 ? "lena.jpg" : argv[1]); // only alloc

    SDL_Surface *postOperationFFT = IMG_Load(argc < 2 ? "lena.jpg" : argv[1]);
    SDL_Surface *postOperationImage = IMG_Load(argc < 2 ? "lena.jpg" : argv[1]);


    /* DO FFT: everything inside the function: */
    //DoFFT(image, imfft);

    SDL_Window *window = SDL_CreateWindow(argv[1],
                                          SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 2 * image->w, image->h, 0);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
    SDL_Texture *texture1 = SDL_CreateTextureFromSurface(renderer, image);
    SDL_Texture *texture2 = SDL_CreateTextureFromSurface(renderer, imfft);


    SDL_Window *postOperationWindow = SDL_CreateWindow("Results",
                                          SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 2 * postOperationImage->w, postOperationImage->h, 0);
    SDL_Renderer *postOperationRenderer = SDL_CreateRenderer(postOperationWindow, -1, 0);
    SDL_Texture *postOperationTextureImage = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationImage);
    SDL_Texture *postOperationTextureFFT = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationFFT);


    datadisplayed_t currentDataDisplayed = datadisplayed_t::DATADISPLAYED_IMAGINARY;
    operationtype_t currentOperationType = operationtype_t::DEFAULT;

    perform_operation(image, imfft, postOperationFFT, postOperationImage, currentOperationType, {}, currentDataDisplayed);

    // redraw the surfaces to textures
    SDL_DestroyTexture(texture1);
    SDL_DestroyTexture(texture2);
    SDL_DestroyTexture(postOperationTextureImage);
    SDL_DestroyTexture(postOperationTextureFFT);

    texture1 = SDL_CreateTextureFromSurface(renderer, image);
    texture2 = SDL_CreateTextureFromSurface(renderer, imfft);
    postOperationTextureImage = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationImage);
    postOperationTextureFFT = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationFFT);


    while (!quit)
    {
        SDL_RenderClear(renderer);
        SDL_RenderClear(postOperationRenderer);
        SDL_WaitEvent(&event);

        switch (event.type)
        {
        case SDL_KEYDOWN:
            handle_input(event, currentOperationType, currentDataDisplayed, quit);
            perform_operation(image, imfft, postOperationFFT, postOperationImage, currentOperationType, {}, currentDataDisplayed);
                
            // redraw the surfaces to textures
            SDL_DestroyTexture(texture1);
            SDL_DestroyTexture(texture2);
            SDL_DestroyTexture(postOperationTextureImage);
            SDL_DestroyTexture(postOperationTextureFFT);

            texture1 = SDL_CreateTextureFromSurface(renderer, image);
            texture2 = SDL_CreateTextureFromSurface(renderer, imfft);
            postOperationTextureImage = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationImage);
            postOperationTextureFFT = SDL_CreateTextureFromSurface(postOperationRenderer, postOperationFFT);
            break;
        case SDL_QUIT:
            quit = true;
            break;
        }

        // render original
        SDL_Rect rect1 = {0, 0, image->w, image->h};
        SDL_RenderCopy(renderer, texture1, NULL, &rect1);
        SDL_Rect rect2 = {image->w, 0, image->w, image->h};
        SDL_RenderCopy(renderer, texture2, NULL, &rect2);
        SDL_RenderPresent(renderer);

        // render post operation effect
        rect1 = {0, 0, postOperationImage->w, postOperationImage->h};
        SDL_RenderCopy(postOperationRenderer, postOperationTextureImage, NULL, &rect1);
        rect2 = {postOperationFFT->w, 0, postOperationFFT->w, postOperationFFT->h};
        SDL_RenderCopy(postOperationRenderer, postOperationTextureFFT, NULL, &rect2);
        SDL_RenderPresent(postOperationRenderer);
    }

    SDL_DestroyTexture(texture1);
    SDL_DestroyTexture(texture2);
    SDL_FreeSurface(image);
    SDL_FreeSurface(imfft);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    SDL_DestroyTexture(postOperationTextureImage);
    SDL_DestroyTexture(postOperationTextureFFT);
    SDL_FreeSurface(postOperationImage);
    SDL_FreeSurface(postOperationFFT);
    SDL_DestroyRenderer(postOperationRenderer);
    SDL_DestroyWindow(postOperationWindow);

    SDL_Quit();
    return 0;
}
