/**
 * Filtro Sobel 
 * 
 * Escrito por: Carlos Peralta Coto
 * 
 * Sin MPI 
**/


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>


/**
 * Esta función aplica el filtro en la matrix img, de proporciones h x w
 * Y devuelve el resultado en imgf de proporciones h-2 x w-2
**/
void filter_y(int h, int w, char img[h][w], int imgf[h-2][w-2])
{
  int maxh = h - 1;
  int maxw = w - 1;
  int result;
  int i;
  int j;
  for (i = 1; i < maxh; i++)
  {
    for (j = 1; j < maxw; j++)
    {
      result = img[i-1][j-1] + 2 * img[i-1][j] + img[i-1][j+1]
               - img[i+1][j-1] - 2 * img[i+1][j] - img[i+1][j+1];
      imgf[i-1][j-1] = result;
    }
  }
}

/**
 * Esta función aplica el filtro en la matrix img, de proporciones h x w
 * Y devuelve el resultado en imgf de proporciones h-2 x w-2
**/
void filter_x(int h, int w, char img[h][w], int imgf[h-2][w-2])
{
  int maxh = h - 1;
  int maxw = w - 1;
  int result;
  int i;
  int j;
  for (i = 1; i < maxh; i++)
  {
    for (j = 1; j < maxw; j++)
    {
      result = img[i-1][j-1] + 2 * img[i][j-1] + img[i+1][j-1]
               - img[i-1][j+1] - 2 * img[i][j+1] - img[i+1][j+1];
      imgf[i-1][j-1] = result;
    }
  }
}


/**
 * Esta función obtiene los valores máximos y mínimos en la matrix imgf, de proporciones h x w
 * Y aplica el mapeo a números de 0 a 255
**/
void map(int h, int w, int imgf[h][w])
{
  int i;
  int j;
  int max = 1443;//imgf[0][0];
  int min = 0;
  
  printf("max = %d, min = %d \n", max, min);
  for (i = 0; i < h; i++)
  {
    for (j = 0; j < w; j++)
    {
      int num = imgf[i][j];
      imgf[i][j] = ((num - min) * 255) / (max - min);
    }
  }
}

/**
 * Función principal, aquí se colocan los valores de height h y width w
 * Abre el archivo, lo filtra, mapea y genera un nuevo archivo con el resultado
**/
int main() 
{

	// estas variables es para medir el tiempo de ejecución del filtrado de la imagen
	clock_t tiempo_inicio, tiempo_final;
	double segundos;
	tiempo_inicio = clock();
	int w, h;

	w = 700;
	h = 500;
	
	FILE *fptr;
	char num[h][w];
	
	//código para el proceso maestro

		
		
		fptr = fopen("./700x500.data","rb"); //Se abre el archivo original
		if(fptr == NULL) 
		{
			printf("Error abriendo archivo");
			return 1;
		}
		//char num[h][w];
		fread(num, sizeof(num), 1, fptr);
		fclose(fptr);
		
	
	int imgf[h-2][w-2];
	int imgf_x[h-2][w-2];
	int imgf_y[h-2][w-2];
	
	filter_x(h,w, num, imgf_x);				//Se filtra en x
	filter_y(h,w, num, imgf_y);				//Se filtra en y

	int i, j;	
	for(i = 0; i < h-2;i++)
	{
		for (j = 0; j < w-2; j++)
		{
			double tmpx = pow(imgf_x[i][j], 2);
			double tmpy = pow(imgf_y[i][j], 2);
			imgf[i][j] = (int)(sqrt(tmpx + tmpy));
		}
	}
	
	map(h-2, w-2, imgf);					//Se mapea la matriz

	fptr = fopen("./filt.data","wb");//Se genera el archivo modificado
	if(fptr == NULL) 
	{
		printf("Error abriendo archivo");
		return 1;
	}
	
	char numf[h-2][w-2];
	for (i = 0; i < h-2; i++)
	{
		for (j = 0; j < w-2; j++)
		{
			numf[i][j] = imgf[i][j];
		}
	}
	fwrite(numf, sizeof(numf), 1, fptr);			//Se guarda el archivo modificado
	fclose(fptr);
	
	tiempo_final = clock();
	segundos = (double)( tiempo_final-tiempo_inicio) / CLOCKS_PER_SEC; // retorna el tiempo en segundos
	printf("El tiempo de ejecución es de: %f \n",segundos); // Imprime el resultado 

	return 0;
}

