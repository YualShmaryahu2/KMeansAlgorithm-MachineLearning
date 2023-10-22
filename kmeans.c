#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Define a linked list structure for individual coordinates.
struct cord
{
    double value;
    struct cord *next;
};

// Define a linked list structure for vectors composed of coordinates.
struct vector
{
    struct vector *next;
    struct cord *cords;
};

// Function prototypes
double **build_table(int num_of_clusters, double **input_data, int dim);
double calc_distance(double *vector1, double *vector2, int dim);
void update_centroid(double **clusters, int* clusters_cnt, int dim, double **new_centroids, int num_of_clusters);
double calc_delta(double **new_centroids, double **prev_centroids, int num_of_clusters, int dim);
void create_output(double **vectors_array, int num_of_clusters, int d);
void kmeans(int K, int iter, double **input_data, int N, int dim);
int check_if_string_is_ok(char* string);
double **createArrayfromInput(int *N ,int *dim_size,struct vector *head,int k,struct cord *a);
double* turn_list_to_array(struct vector *vec, int d);
struct vector free_vector(struct vector *vec);
int number_of_input(struct vector *vectors);
int get_len_from_list(struct vector *vec);
void build_clusters(double **array_of_centroids, double **clusters, int* clusters_cnt, double **input_data, int dim, int N, int num_of_clusters);
void free_matrices(double **array, int num);
void copy(double **array_of_centroids, double **prev_array_of_centroids, int dim, int num_of_clusters);
void freeVectors(struct vector* headVec);

// Function to free a linked list of vectors.
void freeVectors(struct vector* headVec)
{
    struct vector* tmp;
    while (headVec != NULL)
    {
        tmp = headVec;
        headVec = headVec->next;
        free(tmp);
    }
}

// Function to free a linked list of coordinates.
struct vector free_vector(struct vector *vec){
    struct vector *v_next = vec->next;
    struct cord* cordin = vec->cords;
    while (cordin != NULL){
        struct cord* next = cordin->next;
        free (cordin);
        cordin = next;
    }
    free(vec);
    return *v_next;
}

// Function to check if a string contains only numerical characters (0-9).
int check_if_string_is_ok(char* string){
    int i = 0;
    while (string[i] != '\0'){
        int a = string[i];
        if (a < 48 || a > 57){
            return 0;
        }
        i += 1;
    }
    return 1;
}

// Function to free memory allocated for 2D arrays.
void free_matrices(double **array, int num){
    int i;
    i = 0;
    while (i < num){
        free(array[i]);
        i += 1;
    }
    free(array);
}

// Function to count the number of vectors in the linked list.
int number_of_input(struct vector *vectors){
    int i = 0;
    struct vector *ve = vectors;
    while (ve->next != NULL){
        i++;
        ve = ve->next;
    }
    return i;
}

// Function to create a 2D array from the linked list of vectors.
double **createArrayfromInput(int *N, int *dim_size, struct vector *head, int k, struct cord *a){
    struct vector* v;
    double** array_of_vectors;
    double* vector_array;
    int cnt;
    *N = number_of_input(head);
    
    // Check if the number of clusters is valid.
    if (*N <= k || k < 1){
        freeVectors(head);
        free(a);
        printf("Invalid number of clusters!");
        exit(1);
    }
    
    array_of_vectors = (double**)calloc(*N, sizeof(double*));
    cnt = 0;
    
    while (cnt < *N){
        struct vector *s;
        
        if (cnt == 0){
            *dim_size = get_len_from_list(head);
            v = head;
        }
        s = v->next;
        vector_array = turn_list_to_array(v, *dim_size);
        array_of_vectors[cnt] = vector_array;
        cnt += 1;
        free_vector(v);
        v = s;
    }
    return array_of_vectors;
}

// Function to get the dimensionality of a vector.
int get_len_from_list(struct vector *vec){
    int i = 0;
    struct vector *vecs = vec;
    struct cord* cordin = vecs->cords;
    
    while (cordin != NULL){
        i += 1;
        cordin = cordin->next;
    }
    return i;
}

// Function to convert a linked list of coordinates into an array.
double* turn_list_to_array(struct vector *vec, int d){
    int i = 0 ;
    double * val_of_vector = (double*) calloc(d, sizeof(double*));
    struct vector *vecs = vec;
    struct cord* cordin = vecs->cords;
    
    while (i < d){
        val_of_vector[i] = cordin->value;
        cordin = cordin->next;
        i += 1;
    }
    return val_of_vector;
}

// Main K-Means clustering function.
void kmeans(int K, int iter, double **input_data, int N, int dim){
    int i;
    double epsilon = 0.0001;
    int num_of_clusters = K;
    double **array_of_centroids = build_table(num_of_clusters, input_data, dim);
    double **new_centroids = build_table(num_of_clusters, input_data, dim);
    double **temp_array = new_centroids;
    double **temp_array1 = array_of_centroids;
    int iter_cnt = 0;
    double delta = epsilon + 1 ;
    double **prev_array_of_centroids = (double**) malloc(num_of_clusters * sizeof(double*));
    
    for(i = 0; i < num_of_clusters; i++){
        prev_array_of_centroids[i] = (double*) malloc(dim * sizeof(double));
    }
    
    while (iter_cnt < iter && delta >= epsilon){
        double **clusters =(double**) calloc(num_of_clusters, sizeof(double*));
        int *clusters_cnt = (int*) calloc(num_of_clusters, sizeof(int*));
        
        int i;
        iter_cnt ++;
        
        for(i = 0; i < num_of_clusters; i++){
            clusters[i] = (double*) calloc(dim, sizeof(double));
        }
        
        build_clusters(array_of_centroids, clusters, clusters_cnt, input_data, dim, N, num_of_clusters);
        update_centroid(clusters, clusters_cnt, dim, new_centroids, num_of_clusters);
        array_of_centroids = new_centroids;
        
        if (iter_cnt != 1){
            delta = calc_delta(new_centroids, prev_array_of_centroids, num_of_clusters, dim);
        }
        copy(array_of_centroids, prev_array_of_centroids, dim, num_of_clusters);
        free_matrices(clusters, num_of_clusters);
        free(clusters_cnt);
    }
    create_output(array_of_centroids, num_of_clusters, dim);
    free_matrices(temp_array1, num_of_clusters);
    free_matrices(prev_array_of_centroids, num_of_clusters);
    free_matrices(temp_array, num_of_clusters);
}

// Function to copy the contents of one 2D array to another.
void copy(double **array_of_centroids, double **prev_array_of_centroids, int dim, int num_of_clusters){
    int i, j;
    for(i = 0; i < num_of_clusters; i++){
        for(j = 0; j < dim; j++){
            prev_array_of_centroids[i][j] = array_of_centroids[i][j];
        }
    }
}

// Function to build clusters by assigning data points to the nearest centroid.
void build_clusters(double **array_of_centroids, double **clusters, int* clusters_cnt, double **input_data, int dim, int N, int num_of_clusters){
    int i;
    int k;
    for (i = 0; i < N; i++){
        double *vector = input_data[i];
        double min_distance;
        int min_distance_index;
        int j;

        vector = input_data[i];
        min_distance = INFINITY;
        min_distance_index = -1;

        for (j = 0; j < num_of_clusters; j++){
            double *centroid = array_of_centroids[j];
            double curr_distance = calc_distance(vector, centroid, dim);
            if (curr_distance < min_distance){
                min_distance = curr_distance;
                min_distance_index = j;
            }
        }
        clusters_cnt[min_distance_index] += 1;
        for (k = 0; k < dim; k++){
            clusters[min_distance_index][k] += vector[k];
        }
    }
}

// Function to create a 2D array and initialize it with values from another 2D array.
double **build_table(int num_of_clusters, double **input_data, int dim){
    double **array_of_centroids = (double**) malloc(num_of_clusters * sizeof(double*));
    int i, j;
    for(i = 0; i < num_of_clusters; i++){
        array_of_centroids[i] = (double*) malloc(dim * sizeof(double));
        for(j = 0; j < dim; j++){
            array_of_centroids[i][j] = input_data[i][j];
        }
    }
    return array_of_centroids;
}

// Function to calculate the Euclidean distance between two vectors.
double calc_distance(double *vector1, double *vector2, int dim){
    double summ = 0.0;
    int i;
    for(i = 0; i < dim; i++){
        summ += pow((vector1[i] - vector2[i]), 2);
    }
    return sqrt(summ);
}

// Function to update cluster centroids based on assigned data points.
void update_centroid(double **clusters, int* clusters_cnt, int dim, double **new_centroids, int num_of_clusters){
    int i;
    int j;
    double clusters_sum;
    int cluster_num_of_points;
    for (i = 0; i < num_of_clusters; i++){
        for(j = 0; j < dim; j++){
            clusters_sum = clusters[i][j];
            cluster_num_of_points = clusters_cnt[i];
            new_centroids[i][j] =  clusters_sum / cluster_num_of_points;
        }
    }
}

// Function to calculate the maximum change in centroids between iterations.
double calc_delta(double **new_centroids, double **prev_centroids, int num_of_clusters, int dim){
    double delta = 0.0;
    int i;
    for(i = 0; i < num_of_clusters; i++){
        double distance = calc_distance(new_centroids[i], prev_centroids[i], dim);
        if (distance > delta){
            delta = distance;
        }
    }
    return delta;
}

// Function to print the final cluster centroids.
void create_output(double **vectors_array, int num_of_clusters, int dim){
    int i = 0;
    while (i < num_of_clusters){
        int j = 0;
        while (j < dim){
            double num = vectors_array[i][j];
            if (j == dim - 1){
                printf("%.4f\n", num);
            }
            else {
                printf("%.4f,", num);
            }
            j += 1;
        }
        i += 1;
    }
}

int main(int argc, char** argv )
{
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    int check_second_argument, N, k, iter, dim_size;
    char* first_argument;
    char* second_argument;
    double **array_of_vectors;
    double n;
    char c;
    iter = 200;
    check_second_argument = 0;
    
    if (argc > 3 || argc < 2){
        printf("An Error Has Occurred\n");
        return 1;
    }
    else {
        if (argc == 3){
            first_argument = argv[1];
            second_argument = argv[2];
            check_second_argument += 1;
        }
        else {
            first_argument = argv[1];
        }
    }
    
    // Check if the first argument (number of clusters) is a valid integer.
    if (check_if_string_is_ok(first_argument) == 0){
        printf("Invalid number of clusters!");
        return 1;
    }
    
    if (check_second_argument == 1){
        // Check if the second argument (number of iterations) is a valid integer.
        if (check_if_string_is_ok(second_argument) == 0){
            printf("Invalid maximum iteration!");
            return 1;
        }
        iter = atoi(second_argument);
        if (iter <= 1 || iter >= 1000){
            printf("Invalid maximum iteration!");
            return 1;
        }
    }
    
    k = atoi(first_argument);
    
    // Allocate memory for linked lists to read input data.
    head_cord = malloc(sizeof(struct cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;
    
    // Read input data and build linked lists of vectors and coordinates.
    while (scanf("%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
        {
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(struct cord));
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }

        curr_cord->value = n;
        curr_cord->next = malloc(sizeof(struct cord));
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }
    
    // Check if the number of clusters is valid.
    if (check_second_argument == 1){
        if (check_if_string_is_ok(second_argument) == 0){
            printf("Invalid number of clusters!");
            return 1;
        }
        iter = atoi(second_argument);
        if (iter <= 1 || iter >= 1000){
            printf("Invalid maximum iteration!");
            return 1;
        }
        array_of_vectors = createArrayfromInput(&N, &dim_size, head_vec, k, head_cord);
        kmeans(k, iter, array_of_vectors, N, dim_size);
        free_matrices(array_of_vectors, N);
        free(head_cord);
        free(curr_vec);
        return 0;
    }
    else{
        array_of_vectors = createArrayfromInput(&N, &dim_size, head_vec, k, head_cord);
        if (k >= N){
            printf("Invalid number of clusters!");
            return 1;
        }
        kmeans(k, iter, array_of_vectors, N, dim_size);
        free_matrices(array_of_vectors, N);
        free(head_cord);
        free(curr_vec);
        return 0;
    }
}
