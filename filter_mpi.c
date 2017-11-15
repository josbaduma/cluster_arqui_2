/**
 * Filtro Sobel 
 * 
 * Escrito por: Carlos Peralta Coto
 * 
 * Implementa MPI para la comunicación serial
**/

#include <mpi.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/*MPI_Send(
    void* data,
    int count,
    MPI_Datatype datatype,
    int destination,
    int tag,
    MPI_Comm communicator)
MPI_Recv(
    void* data,
    int count,
    MPI_Datatype datatype,
    int source,
    int tag,
    MPI_Comm communicator,
    MPI_Status* status)
    */
/**
 * Esta función aplica el filtro en la matrix img, de proporciones h x w
 * Y devuelve el resultado en imgf de proporciones h-2 x w-2
**/
void filter_y(int h, int w, char img[h][w], int imgf[h-2][w-2], int inicio, int final)
{
  int maxh = h - 1;
  int maxw = w - 1;
  int result;
  int i;
  int j;
  if(inicio == 0){inicio = 1;}
  if(final == h){final = h-1;}
  for (i = inicio; i < final; i++)
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
void filter_x(int h, int w, char img[h][w], int imgf[h-2][w-2], int inicio, int final)
{
  int maxh = h - 1;
  int maxw = w - 1;
  int result;
  int i;
  int j;
  if(inicio == 0){inicio = 1;}
  if(final == h){final = h-1;}
  for (i = inicio; i < final; i++)
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
  /*for (i = 0; i < h; i++)
  {
    for (j = 0; j < w; j++)
    {
      int num = imgf[i][j];
      max = (num > max)? num:max;
      min = (num < min)? num:min;
    }
  }*/
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


	int matriz[10][10];
	int fila[10];
	for (int i=0;i<10;i++){
		for(int j=0;j<10;j++){
			matriz[i][j]=i*j;
		}
	}

	int world_size, rank, name_len, ierr, w, h, root_process;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	ierr = MPI_Init(NULL, NULL);// iniciocializa el ambiente MPI
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &world_size);// Guarda la cantidad de procesos
   	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);// Guarda el número del rank
	MPI_Get_processor_name(processor_name, &name_len);// Guarda el nombre del procesador
	MPI_Status Stat;

	root_process = 0;
	w = 646;
	h = 443;
	
	FILE *fptr;
	char num[h][w];
	int imagenfiltrada[h][w-2];
	//código para el proceso maestro
	if(rank == root_process)
	{
		
		
		fptr = fopen("./646x443.data","rb"); //Se abre el archivo original
		if(fptr == NULL) 
		{
			printf("Error abriendo archivo");
			return 1;
		}
		//char num[h][w];
		fread(num, sizeof(num), 1, fptr);
		fclose(fptr);

		int filapornodo=(int) (h/(world_size-1));

		for(int i=1;i<world_size;i++){
			int inicio=(i-1)*filapornodo;
			int final=i*filapornodo;
			MPI_Send(&inicio,1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&final,1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&num[0][0],h*w, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		}

		for(int i=1;i<world_size;i++){
			int inicio=(i-1)*filapornodo;
			int final=i*filapornodo;
			
			int matriz[final-inicio][w-2];
			//int imagenfiltrada[h][w-2];
			MPI_Recv(&matriz, (final-inicio)*(w-2), MPI_INT, i, 1, MPI_COMM_WORLD,&Stat);
			for(int k = 0; k < (final-inicio);k++)
			{
				for (int x = 0; x < w-2; x++)
				{
					
					imagenfiltrada[inicio+k][x]=matriz[k][x];
				}
				
			}
		}
		printf("\nSe envió %d \n",world_size);
		
		fptr = fopen("./filt.data","wb");//Se genera el archivo modificado
		if(fptr == NULL) 
		{
			printf("Error abriendo archivo");
			return 1;
		}
		
		char numf[h-2][w-2];
		for (int i = 0; i < h-2; i++)
		{
			for (int j = 0; j < w-2; j++)
			{
				numf[i][j] = imagenfiltrada[i][j];
			}
		}
		fwrite(numf, sizeof(numf), 1, fptr);//Se guarda el archivo modificado
		fclose(fptr);
		//variable que toma el timpo al finalizar la ejecución

		tiempo_final = clock();
		segundos = (double)( tiempo_final-tiempo_inicio) / CLOCKS_PER_SEC; // retorna el tiempo en segundos
		printf("El tiempo de ejecución es de: %f \n",segundos); // Imprime el resultado 

	}
	else {
		
		char fila[h][w];
		int inicio=0;
		int final=0;
			MPI_Recv(&inicio, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&Stat);
			MPI_Recv(&final, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&Stat);

			MPI_Recv(&num, h*w, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&Stat);
//			printf("Soy proceso %d del procesador: %s y analizo desde %d hasta %d \n",rank, processor_name,inicio,final);

			int imgf[h-2][w-2];
			int imgf_x[h-2][w-2];
			int imgf_y[h-2][w-2];
			
			filter_x(h,w, num, imgf_x, inicio, final);				//Se filtra en x
			filter_y(h,w, num, imgf_y, inicio, final);			    //Se filtra en y

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

			int altura = final - inicio;
			if (inicio == 0) altura = altura - 1;
			if (final == h-1) altura = altura -1;
			int imagen[altura][w-2];
	
			for(i = 0; i < altura;i++)
			{
				for(j = 0; j < w-2; j++)
					//imagen[i][j] = imgf[inicio+i][j];
					imagen[i][j] = (imgf[inicio+i][j] * 255) / (1443);
					//printf("%d,",num[i][j]);
			}
			MPI_Send(&imagen[0][0],altura*(w-2), MPI_INT, 0, 1, MPI_COMM_WORLD);
}
	MPI_Finalize(); // finalaliza el ambiente MPI.
	return 0;

}

