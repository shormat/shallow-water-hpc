#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

struct parameters
{
    double dx, dy, dt, max_t;
    double g, gamma;
    int source_type;
    int sampling_rate;
    char input_h_filename[256];
    char output_eta_filename[256];
    char output_u_filename[256];
    char output_v_filename[256];
};

struct data
{
    int nx, ny;
    int offsetx, offsety;
    bool has_ghost;
    int ld;
    double dx, dy;
    double *values;
    MPI_Datatype row_t, col_t;
};

typedef struct process_t
{
    MPI_Comm comm;
    int rank;
    int coords[2];          // [cy, cx] (Y, X)
    int neighbors[4];       // UP, DOWN, LEFT, RIGHT
    int offset_x, offset_y; // global interior start indices
    int local_nx, local_ny; // interior sizes
} Process;

enum directions
{
    UP,
    DOWN,
    LEFT,
    RIGHT
};

#define GET(D,i,j) ((D)->values[((j) + (D)->has_ghost) * (D)->ld + ((i) + (D)->has_ghost)])
#define SET(D,i,j,val) (GET((D),(i),(j)) = (val))

int read_parameters(struct parameters *param, const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        printf("Error: Could not open parameter file '%s'\n", filename);
        return 1;
    }
    int ok = 1;
    if (ok)
        ok = (fscanf(fp, "%lf", &param->dx) == 1);
    if (ok)
        ok = (fscanf(fp, "%lf", &param->dy) == 1);
    if (ok)
        ok = (fscanf(fp, "%lf", &param->dt) == 1);
    if (ok)
        ok = (fscanf(fp, "%lf", &param->max_t) == 1);
    if (ok)
        ok = (fscanf(fp, "%lf", &param->g) == 1);
    if (ok)
        ok = (fscanf(fp, "%lf", &param->gamma) == 1);
    if (ok)
        ok = (fscanf(fp, "%d", &param->source_type) == 1);
    if (ok)
        ok = (fscanf(fp, "%d", &param->sampling_rate) == 1);
    if (ok)
        ok = (fscanf(fp, "%255s", param->input_h_filename) == 1);
    if (ok)
        ok = (fscanf(fp, "%255s", param->output_eta_filename) == 1);
    if (ok)
        ok = (fscanf(fp, "%255s", param->output_u_filename) == 1);
    if (ok)
        ok = (fscanf(fp, "%255s", param->output_v_filename) == 1);
    fclose(fp);
    if (!ok)
    {
        printf("Error: Could not read one or more parameters in '%s'\n", filename);
        return 1;
    }
    return 0;
}

void print_parameters(const struct parameters *param)
{
    printf("Parameters:\n");
    printf(" - grid spacing (dx, dy): %g m, %g m\n", param->dx, param->dy);
    printf(" - time step (dt): %g s\n", param->dt);
    printf(" - maximum time (max_t): %g s\n", param->max_t);
    printf(" - gravitational acceleration (g): %g m/s^2\n", param->g);
    printf(" - dissipation coefficient (gamma): %g 1/s\n", param->gamma);
    printf(" - source type: %d\n", param->source_type);
    printf(" - sampling rate: %d\n", param->sampling_rate);
    printf(" - input bathymetry (h) file: '%s'\n", param->input_h_filename);
    printf(" - output elevation (eta) file: '%s'\n", param->output_eta_filename);
    printf(" - output velocity (u, v) files: '%s', '%s'\n",
           param->output_u_filename, param->output_v_filename);
}

int read_data(struct data *data, const char *filename)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
    {
        printf("Error: Could not open input data file '%s'\n", filename);
        return 1;
    }
    int ok = 1;
    if (ok)
        ok = (fread(&data->nx, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->ny, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->dx, sizeof(double), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->dy, sizeof(double), 1, fp) == 1);
    if (ok)
    {
        int N = data->nx * data->ny;
        if (N <= 0)
        {
            printf("Error: Invalid number of data points %d\n", N);
            ok = 0;
        }
        else
        {
            data->values = (double *)malloc(N * sizeof(double));
            if (!data->values)
            {
                printf("Error: Could not allocate data (%d doubles)\n", N);
                ok = 0;
            }
            else
            {
                ok = (fread(data->values, sizeof(double), N, fp) == N);
            }
        }
    }
    fclose(fp);
    if (!ok)
    {
        printf("Error reading input data file '%s'\n", filename);
        return 1;
    }

    
    data->has_ghost = false;
    data->ld = data->nx;
    data->offsetx = 0;
    data->offsety = 0;
    data->row_t = MPI_DATATYPE_NULL;
    data->col_t = MPI_DATATYPE_NULL;
    
    return 0;
}

int write_data(const struct data *data, const char *filename, int step)
{
    char out[512];
    if (step < 0)
        sprintf(out, "%s.dat", filename);
    else
        sprintf(out, "%s_%d.dat", filename, step);
    FILE *fp = fopen(out, "wb");
    if (!fp)
    {
        printf("Error: Could not open output data file '%s'\n", out);
        return 1;
    }
    int ok = 1;
    if (ok)
        ok = (fwrite(&data->nx, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fwrite(&data->ny, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fwrite(&data->dx, sizeof(double), 1, fp) == 1);
    if (ok)
        ok = (fwrite(&data->dy, sizeof(double), 1, fp) == 1);
    int N = data->nx * data->ny;
    if (ok)
        ok = (fwrite(data->values, sizeof(double), N, fp) == N);
    fclose(fp);
    if (!ok)
    {
        printf("Error writing data file '%s'\n", out);
        return 1;
    }
    return 0;
}

int write_data_vtk(const struct data *data, const char *name,
                   const char *filename, int step, int rank)
{
  char out[512];
  if(step < 0)
    sprintf(out, "%s_rank%d.vti", filename, rank);
  else
    sprintf(out, "%s_rank%d_%d.vti", filename, rank, step);

  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output VTK file '%s'\n", out);
    return 1;
  }

  uint64_t num_points = data->nx * data->ny;
  uint64_t num_bytes = num_points * sizeof(double);

  fprintf(fp, "<?xml version=\"1.0\"?>\n"
              "<VTKFile"
              " type=\"ImageData\""
              " version=\"1.0\""
              " byte_order=\"LittleEndian\""
              " header_type=\"UInt64\""
              ">\n"
              "  <ImageData"
              " WholeExtent=\"0 %d 0 %d 0 0\""
              " Spacing=\"%lf %lf 0.0\""
              " Origin=\"%lf %lf 0\""
              ">\n"
              "    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
              "      <PointData Scalars=\"scalar_data\">\n"
              "        <DataArray"
              " type=\"Float64\""
              " Name=\"%s\""
              " format=\"appended\""
              " offset=\"0\""
              ">\n"
              "        </DataArray>\n"
              "      </PointData>\n"
              "    </Piece>\n"
              "  </ImageData>\n"
              "  <AppendedData encoding=\"raw\">\n_",
              data->nx - 1, data->ny - 1,
              data->dx, data->dy,
              data->offsetx * data->dx, data->offsety * data->dy,
              data->nx - 1, data->ny - 1,
              name);

  fwrite(&num_bytes, sizeof(uint64_t), 1, fp);
  fwrite(data->values, sizeof(double), num_points, fp);

  fprintf(fp, "  </AppendedData>\n"
              "</VTKFile>\n");

  fclose(fp);

  return 0;
}

int write_manifest_vtk(const char *filename, double dt, int nt,
                       int sampling_rate, int numranks)
{
  char out[512];
  sprintf(out, "%s.pvd", filename);

  FILE *fp = fopen(out, "wb");

  if(!fp) {
    printf("Error: Could not open output VTK manifest file '%s'\n", out);
    return 1;
  }

  fprintf(fp, "<VTKFile"
              " type=\"Collection\""
              " version=\"0.1\""
              " byte_order=\"LittleEndian\">\n"
              "  <Collection>\n");

  for(int n = 0; n < nt; n++) {
    if(sampling_rate && !(n % sampling_rate)) {
      double t = n * dt;
      for (int rank = 0; rank < numranks; rank++) {
        fprintf(fp, "    <DataSet timestep=\"%g\" part=\"%d\" file='%s_rank%d_%d.vti'/>\n",
                t, rank, filename, rank, n);
      }
    }
  }

  fprintf(fp, "  </Collection>\n"
              "</VTKFile>\n");

  fclose(fp);

  return 0;
}


int init_data(struct data *data, int nx, int ny, int offsetx, int offsety,
              double dx, double dy, double val)
{
    data->nx = nx;
    data->ny = ny;
    data->offsetx = offsetx;
    data->offsety = offsety;
    data->dx = dx;
    data->dy = dy;
    data->has_ghost = true;
    data->ld = nx + 2 * data->has_ghost;

    size_t total_size = data->ld * (ny + 2 * data->has_ghost);
    data->values = (double *)malloc(total_size * sizeof(double));
    if (!data->values)
    {
        printf("Error: Could not allocate data\n");
        return 1;
    }
    #pragma omp parallel for
    for(size_t i = 0; i < total_size; i++)
        data->values[i] = val;

    MPI_Type_contiguous(nx, MPI_DOUBLE, &data->row_t);  MPI_Type_commit(&data->row_t);
    MPI_Type_vector(ny, 1, data->ld, MPI_DOUBLE, &data->col_t); MPI_Type_commit(&data->col_t);
    return 0;
}

void free_data(struct data *d)
{   
    if (d->values) free(d->values);
    if (d->row_t != MPI_DATATYPE_NULL) MPI_Type_free(&d->row_t);
    if (d->col_t != MPI_DATATYPE_NULL) MPI_Type_free(&d->col_t);
    d->values = NULL;
    d->row_t = d->col_t = MPI_DATATYPE_NULL;
}

double interpolate_data(const struct data *src_no_ghost, double x, double y)
{
  const struct data *data = src_no_ghost; // clarity
  double gx = x / data->dx;
  double gy = y / data->dy;

  int ix = (int)floor(gx);
  int iy = (int)floor(gy);

  if (ix < 0) ix = 0; else if (ix >= data->nx - 1) ix = data->nx - 2;
  if (iy < 0) iy = 0; else if (iy >= data->ny - 1) iy = data->ny - 2;

  // (capped at 1 so points on the far edge take the edge value)
  double tx = fmin(gx - ix, 1.);
  double ty = fmin(gy - iy, 1.);

  double v00 = data->values[iy * data->nx + ix];
  double v10 = data->values[iy * data->nx + (ix+1)];
  double v01 = data->values[(iy+1) * data->nx + ix];
  double v11 = data->values[(iy+1) * data->nx + (ix+1)];

  return v00*(1-tx)*(1-ty) + v10*tx*(1-ty) + v01*(1-tx)*ty + v11*tx*ty;
}

static inline int block_start(int N, int P, int i) { return (i * N) / P; }
static inline int block_size(int N, int P, int i) { return ((i + 1) * N) / P - (i * N) / P; }

Process *init_process(MPI_Comm cart_comm, int Nx, int Ny, int Px, int Py)
{
    Process *p = (Process*) malloc(sizeof(Process));
    if(!p)
        return NULL;

    p->comm = cart_comm;
    MPI_Comm_rank(cart_comm, &p->rank);
    MPI_Cart_coords(cart_comm, p->rank, 2, p->coords);

    // Axis 0 = Y (rows), axis 1 = X (columns)
    MPI_Cart_shift(cart_comm, 0, 1, &p->neighbors[UP], &p->neighbors[DOWN]);
    MPI_Cart_shift(cart_comm, 1, 1, &p->neighbors[LEFT], &p->neighbors[RIGHT]);

    int coord_y = p->coords[0];
    int coord_x = p->coords[1];

    p->offset_y = block_start(Ny, Py, coord_y);
    p->offset_x = block_start(Nx, Px, coord_x);
    p->local_ny = block_size(Ny, Py, coord_y);
    p->local_nx = block_size(Nx, Px, coord_x);

    return p;
}

void free_process(Process *proc)
{
    if (proc) 
        free(proc);
}

/*
 * The explicit scheme only reads one neighbour per field and direction:
 *   eta update reads u(i+1, j) and v(i, j+1)   -> ghost column nx of u, row ny of v
 *   u update reads eta(i-1, j)                 -> ghost column -1 of eta
 *   v update reads eta(i, j-1)                 -> ghost row -1 of eta
 * so each ghost layer is filled from a single neighbour.
 */
static inline void fill_ghost_col_right(const Process *p, struct data *d)
{
    MPI_Sendrecv(&GET(d, 0, 0), 1, d->col_t, p->neighbors[LEFT], 10,
                 &GET(d, d->nx, 0), 1, d->col_t, p->neighbors[RIGHT], 10,
                 p->comm, MPI_STATUS_IGNORE);
}

static inline void fill_ghost_col_left(const Process *p, struct data *d)
{
    MPI_Sendrecv(&GET(d, d->nx - 1, 0), 1, d->col_t, p->neighbors[RIGHT], 11,
                 &GET(d, -1, 0), 1, d->col_t, p->neighbors[LEFT], 11,
                 p->comm, MPI_STATUS_IGNORE);
}

static inline void fill_ghost_row_top(const Process *p, struct data *d)
{
    MPI_Sendrecv(&GET(d, 0, 0), 1, d->row_t, p->neighbors[UP], 20,
                 &GET(d, 0, d->ny), 1, d->row_t, p->neighbors[DOWN], 20,
                 p->comm, MPI_STATUS_IGNORE);
}

static inline void fill_ghost_row_bottom(const Process *p, struct data *d)
{
    MPI_Sendrecv(&GET(d, 0, d->ny - 1), 1, d->row_t, p->neighbors[DOWN], 21,
                 &GET(d, 0, -1), 1, d->row_t, p->neighbors[UP], 21,
                 p->comm, MPI_STATUS_IGNORE);
}

void best_dim(int num_processes, int Nx, int Ny, int dims[2])
{
    double Px_star = sqrt((double)num_processes * Nx / Ny);
    long best_Px = 1;
    long best_Py = num_processes;
    long best_cost = (long)Ny * best_Px + (long)Nx * best_Py;
    for (long px = 1; px * px <= num_processes; ++px)
    {
        if (num_processes % px)
            continue;
        long py = num_processes / px;
        long cost1 = (long)Ny * px + (long)Nx * py;
        if (cost1 < best_cost)
        {
            best_Px = px;
            best_Py = py;
            best_cost = cost1;
        }
        long cost2 = (long)Ny * py + (long)Nx * px;
        if (cost2 < best_cost)
        {
            best_Px = py;
            best_Py = px;
            best_cost = cost2;
        }
    }
    dims[0] = (int)best_Px;
    dims[1] = (int)best_Py;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        printf("Usage: %s parameter_file\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int world_size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank)
    {
        printf(" - MPI processes: %d\n", world_size);
    }

    struct parameters param;
    if (read_parameters(&param, argv[1]))
    {
        MPI_Finalize();
        return 1;
    }

    if(!rank)
        print_parameters(&param);

    struct data h;
    if (read_data(&h, param.input_h_filename))
    { MPI_Finalize(); return 1; }

    // infer size of domain from input bathymetric data
    double hx = h.nx * h.dx;
    double hy = h.ny * h.dy;
    int nx = floor(hx / param.dx);
    int ny = floor(hy / param.dy);
    if (nx <= 0)
        nx = 1;
    if (ny <= 0)
        ny = 1;
    int nt = floor(param.max_t / param.dt);

    if (!rank)
    {
        printf(" - grid size: %g m x %g m (%d x %d = %d grid points)\n",
               hx, hy, nx, ny, nx * ny);
        printf(" - number of time steps: %d\n", nt);
    }

    
    int dims_best[2]; // {Px, Py}
    int periods[2] = {0, 0};
    best_dim(world_size, nx, ny, dims_best);
    int Px = dims_best[0], Py = dims_best[1];

    int dims_cart[2] = {Py, Px}; // [Y, X]
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims_cart, periods, 0, &cart_comm);

    int cart_rank = -1;
    MPI_Comm_rank(cart_comm, &cart_rank);

    Process *proc = init_process(cart_comm, nx, ny, Px, Py);
    if (!proc)
    {
        printf("Error: Could not initialize process %d\n", cart_rank);
        MPI_Finalize();
        return 1;
    }

    
    struct data eta, u, v, h_interp;
    if (init_data(&eta, proc->local_nx, proc->local_ny,
                  proc->offset_x, proc->offset_y, param.dx, param.dy, 0.0))
    {
        MPI_Finalize();
        return 1;
    }
    if (init_data(&u, proc->local_nx, proc->local_ny,
                  proc->offset_x, proc->offset_y, param.dx, param.dy, 0.0))
    {
        MPI_Finalize();
        return 1;
    }
    if (init_data(&v, proc->local_nx, proc->local_ny,
                  proc->offset_x, proc->offset_y, param.dx, param.dy, 0.0))
    {
        MPI_Finalize();
        return 1;
    }
    if (init_data(&h_interp, proc->local_nx, proc->local_ny,
                  proc->offset_x, proc->offset_y, param.dx, param.dy, 0.0))
    {
        MPI_Finalize();
        return 1;
    }

    
    // Interpolate the right ghost column and top ghost row too (read by the
    // eta update): same value as the neighbour would send, and the true depth
    // on the global boundary instead of 0
    #pragma omp parallel for collapse(2)
    for (int j = 0; j <= proc->local_ny; ++j)
    {
        for (int i = 0; i <= proc->local_nx; ++i)
        {
            double x = (proc->offset_x + i) * param.dx;
            double y = (proc->offset_y + j) * param.dy;
            double val = interpolate_data(&h, x, y);
            SET(&h_interp, i, j, val);
        }
    }

    double start = GET_TIME();

    // Constants
    const double A = 5;
    const double f = 1. / 20.;

    // Constants for update
    const double c1 = param.dt * param.g;
    const double c2 = param.dt * param.gamma;
    const double one_minus_c2 = 1. - c2;
    const double c1_over_dx = c1 / param.dx;
    const double c1_over_dy = c1 / param.dy;

    const int AT_LEFT = (proc->neighbors[LEFT] == MPI_PROC_NULL);
    const int AT_RIGHT = (proc->neighbors[RIGHT] == MPI_PROC_NULL);
    const int AT_BOTTOM = (proc->neighbors[UP] == MPI_PROC_NULL);
    const int AT_TOP = (proc->neighbors[DOWN] == MPI_PROC_NULL);

    double *temp_write = (double *)malloc(eta.nx * eta.ny * sizeof(double));
    if(!temp_write)
    {
        printf("Error: Could not allocate temporary write buffer\n");
        free_data(&eta);
        free_data(&u);
        free_data(&v);
        free_data(&h_interp);
        free_process(proc);
        MPI_Finalize();
        return 1;
    }

    #pragma omp parallel
    {
        for (int n = 0; n < nt; n++)
        {

        #pragma omp single
        {
            if (n && (n % (nt / 10)) == 0)
            {
                double time_sofar = GET_TIME() - start;
                double eta = (nt - n) * time_sofar / n;
                printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
                fflush(stdout);
            }
        }

        
        if (param.sampling_rate && !(n % param.sampling_rate))
        {
            #pragma omp for collapse(2)
            for (int j = 0; j < eta.ny; ++j)
                for (int i = 0; i < eta.nx; ++i)
                    temp_write[(size_t)j * eta.nx + i] = GET(&eta, i, j);

            #pragma omp single
            {
                struct data view = eta; // shallow copy for writer
                view.has_ghost = false;
                view.ld = eta.nx;   // contiguous view
                view.values = temp_write; // packed interior

                write_data_vtk(&view, "water elevation", param.output_eta_filename, n, rank);
                // write_data_vtk(&u, "x velocity", param.output_u_filename, n);
                // write_data_vtk(&v, "y velocity", param.output_v_filename, n);
            }
        }

        // impose boundary conditions
        double t = n * param.dt;
        double A_sin = A * sin(2 * M_PI * f * t);
        if (param.source_type == 1)
        {
            /* Example: no-through-flow on left/right/bottom; top forcing on v */
            if (AT_LEFT)
                #pragma omp for 
                for (int j = 0; j < u.ny; ++j)
                    SET(&u, 0, j, 0.0);
            if (AT_RIGHT)
                #pragma omp for 
                for (int j = 0; j < u.ny; ++j)
                    SET(&u, u.nx, j, 0.0);
            if (AT_BOTTOM)
                #pragma omp for 
                for (int i = 0; i < v.nx; ++i)
                    SET(&v, i, 0, 0.0);
            if (AT_TOP)
                #pragma omp for 
                for (int i = 0; i < v.nx; ++i)
                    SET(&v, i, v.ny, A_sin);
        }
        else if (param.source_type == 2)
        {
            int gx = nx / 2, gy = ny / 2;
            int li = gx - proc->offset_x;
            int lj = gy - proc->offset_y;
            if (0 <= li && li < eta.nx && 0 <= lj && lj < eta.ny)
            {
                SET(&eta, li, lj, A_sin);
            }
        }
        else
        {
            if (!rank)
                printf("Error: Unknown source type %d\n", param.source_type);
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        #pragma omp barrier
        #pragma omp master
        {
            fill_ghost_col_right(proc, &u);
            fill_ghost_row_top(proc, &v);
        }
        #pragma omp barrier


        // update eta
        #pragma omp for collapse(2)
        for (int j = 0; j < eta.ny; j++)
        {
            for (int i = 0; i < eta.nx; i++)
            {
                double h_u = (GET(&h_interp, i, j) + GET(&h_interp, i + 1, j)) / 2.0;
                double h_v = (GET(&h_interp, i, j) + GET(&h_interp, i, j + 1)) / 2.0;

                double c1_u = param.dt * h_u;
                double c1_v = param.dt * h_v;

                double eta_ij = GET(&eta, i, j) -
                                c1_u / param.dx * (GET(&u, i + 1, j) - GET(&u, i, j)) -
                                c1_v / param.dy * (GET(&v, i, j + 1) - GET(&v, i, j));

                SET(&eta, i, j, eta_ij);
            }
        }

        #pragma omp master
        {
            fill_ghost_col_left(proc, &eta);
            fill_ghost_row_bottom(proc, &eta);
        }
        #pragma omp barrier

        // update u and v
        #pragma omp for collapse(2) nowait
        for (int j = 0; j < u.ny; ++j)
        {
            for (int i = 0; i < u.nx; ++i)
            {
                double u_ij = one_minus_c2 * GET(&u, i, j) - c1_over_dx * (GET(&eta, i, j) - GET(&eta, i - 1, j));
                SET(&u, i, j, u_ij);
            }
        }
        
        #pragma omp for collapse(2)
        for (int j = 0; j < v.ny; ++j)
        {
            for (int i = 0; i < v.nx; ++i)
            {
                double v_ij = one_minus_c2 * GET(&v, i, j) - c1_over_dy * (GET(&eta, i, j) - GET(&eta, i, j - 1));
                SET(&v, i, j, v_ij);
            }
        }


    }
    }

    if(!rank)
        write_manifest_vtk(param.output_eta_filename, param.dt, nt,
                           param.sampling_rate, world_size);
    // write_manifest_vtk(param.output_u_filename, param.dt, nt,
    //                    param.sampling_rate, world_size);
    // write_manifest_vtk(param.output_v_filename, param.dt, nt,
    //                    param.sampling_rate, world_size);

    double time = GET_TIME() - start;
    if(!rank)
        printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
               1e-6 * (double)eta.nx * (double)eta.ny * (double)nt / time);

    free_data(&h_interp);
    free_data(&eta);
    free_data(&u);
    free_data(&v);
    free_data(&h);
    free_process(proc);
    free(temp_write);
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}