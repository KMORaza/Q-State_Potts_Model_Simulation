#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <glib.h>

#define MAX_Q 20
#define MAX_DIM 3
#define MAX_SIZE 256
#define HIGHLIGHT_DURATION_MS 200 
#define CELL_SIZE 10 

typedef struct {
    int q;                  // Number of states (Q)
    int dim;                // Dimension (2 or 3)
    int size[MAX_DIM];      // Lattice size in each dimension
    int geometry;           // 0=square, 1=triangular, 2=hexagonal
    int boundary;           // 0=periodic, 1=open
    double temperature;     // Simulation temperature
    int algorithm;          // 0=Metropolis, 1=Heat Bath, 2=Swendsen-Wang
    int initialized;        // Whether lattice is initialized
    int steps;              // Number of steps performed
    int *lattice;           // The lattice itself
    double energy;          // Current energy
    double magnetization;   // Current magnetization
    int *changed_sites;     // Array to track changed sites
    int changed_count;      // Number of changed sites
    gint64 *change_times;   // Timestamps for highlighting changed sites
} PottsModel;

// Cluster information for Swendsen-Wang
typedef struct {
    int *labels;            // Cluster labels for each site
    int *parent;            // Union-Find parent array
    int *size;              // Size of each cluster
    int total_clusters;     // Total number of clusters
} ClusterInfo;

typedef struct {
    GtkWidget *window;
    GtkWidget *drawing_area;
    GtkWidget *speed_scale;
    GtkWidget *start_button;
    GtkWidget *stop_button;
    GtkWidget *pause_button;
    GtkWidget *restart_button;
    GtkWidget *q_scale;
    GtkWidget *size_scale;
    GtkWidget *temp_scale;
    GtkWidget *combo_algorithm;
    GtkWidget *label_energy;
    GtkWidget *label_magnetization;
    GtkWidget *label_steps;
    PottsModel model;
    ClusterInfo clusters;
    bool running;
    bool paused;
    double speed;           
    guint timeout_id;      
} AppData;

GdkRGBA colors[MAX_Q];

void init_colors();
int find_root(int *parent, int x);
void union_clusters(int *parent, int *size, int x, int y);
void potts_init(PottsModel *model, ClusterInfo *clusters);
void potts_calculate_energy(PottsModel *model);
void potts_calculate_magnetization(PottsModel *model);
void potts_metropolis_step(PottsModel *model);
void potts_heat_bath_step(PottsModel *model);
void potts_swendsen_wang_step(PottsModel *model, ClusterInfo *clusters);
void potts_step(PottsModel *model, ClusterInfo *clusters);
void draw_callback(GtkDrawingArea *area, cairo_t *cr, int width, int height, gpointer user_data);
void start_simulation(GtkButton *button, gpointer user_data);
void stop_simulation(GtkButton *button, gpointer user_data);
void pause_simulation(GtkButton *button, gpointer user_data);
void restart_simulation(GtkButton *button, gpointer user_data);
void speed_changed(GtkRange *range, gpointer user_data);
void activate(GtkApplication *app, gpointer user_data);
void init_colors() {
    for (int i = 0; i < MAX_Q; i++) {
        double hue = (double)i / MAX_Q;
        colors[i].red = 0.5 + 0.5 * sin(2 * M_PI * hue);
        colors[i].green = 0.5 + 0.5 * sin(2 * M_PI * hue + 2 * M_PI / 3);
        colors[i].blue = 0.5 + 0.5 * sin(2 * M_PI * hue + 4 * M_PI / 3);
        colors[i].alpha = 1.0;
    }
}

/// Union-Find operations for Swendsen-Wang
int find_root(int *parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]]; 
        x = parent[x];
    }
    return x;
}

void union_clusters(int *parent, int *size, int x, int y) {
    int root_x = find_root(parent, x);
    int root_y = find_root(parent, y);
    if (root_x == root_y) return;
    /// Union by size
    if (size[root_x] < size[root_y]) {
        parent[root_x] = root_y;
        size[root_y] += size[root_x];
    } else {
        parent[root_y] = root_x;
        size[root_x] += size[root_y];
    }
}

void potts_init(PottsModel *model, ClusterInfo *clusters) {
    if (model->initialized) {
        free(model->lattice);
        free(model->changed_sites);
        free(model->change_times);
        if (clusters->labels) {
            free(clusters->labels);
            free(clusters->parent);
            free(clusters->size);
        }
    }
    int total_size = 1;
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
    }
    model->lattice = malloc(total_size * sizeof(int));
    model->changed_sites = malloc(total_size * sizeof(int));
    model->change_times = malloc(total_size * sizeof(gint64));
    clusters->labels = malloc(total_size * sizeof(int));
    clusters->parent = malloc(total_size * sizeof(int));
    clusters->size = malloc(total_size * sizeof(int));
    model->initialized = 1;
    model->steps = 0;
    model->changed_count = 0;
    for (int i = 0; i < total_size; i++) {
        model->lattice[i] = rand() % model->q;
        model->change_times[i] = 0;
    }
    for (int i = 0; i < total_size; i++) {
        clusters->parent[i] = i;
        clusters->size[i] = 1;
    }
    
    /// Calculate initial energy and magnetization
    potts_calculate_energy(model);
    potts_calculate_magnetization(model);
}

/// Calculate the energy of the current configuration
void potts_calculate_energy(PottsModel *model) {
    double energy = 0.0;
    int total_size = 1;
    int strides[MAX_DIM] = {1};
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
        if (d > 0) {
            strides[d] = strides[d-1] * model->size[d-1];
        }
    }
    for (int i = 0; i < total_size; i++) {
        int current_state = model->lattice[i];
        for (int d = 0; d < model->dim; d++) {
            int neighbor_pos = i + strides[d];
            if ((i / strides[d]) % model->size[d] == model->size[d] - 1) {
                if (model->boundary == 0) { // Periodic
                    neighbor_pos -= model->size[d] * strides[d];
                } else { // Open
                    continue;
                }
            }
            if (neighbor_pos < total_size) {
                energy += (current_state == model->lattice[neighbor_pos]) ? -1.0 : 0.0;
            }
            // Negative neighbor (for non-directional interactions)
            int neighbor_neg = i - strides[d];
            // Handle boundary conditions
            if ((i / strides[d]) % model->size[d] == 0) {
                if (model->boundary == 0) { // Periodic
                    neighbor_neg += model->size[d] * strides[d];
                } else { // Open
                    continue;
                }
            }
            if (neighbor_neg >= 0) {
                energy += (current_state == model->lattice[neighbor_neg]) ? -1.0 : 0.0;
            }
        }
    }
    model->energy = energy;
}

/// Calculate the magnetization of the current configuration
void potts_calculate_magnetization(PottsModel *model) {
    int total_size = 1;
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
    }
    // Q*max(n_i) - N)/(N*(Q-1))
    int *state_counts = calloc(model->q, sizeof(int));
    for (int i = 0; i < total_size; i++) {
        state_counts[model->lattice[i]]++;
    }
    int max_count = 0;
    for (int i = 0; i < model->q; i++) {
        if (state_counts[i] > max_count) {
            max_count = state_counts[i];
        }
    }
    model->magnetization = (model->q * max_count - total_size) / (double)(total_size * (model->q - 1));
    free(state_counts);
}

/// Metropolis algorithm step
void potts_metropolis_step(PottsModel *model) {
    int total_size = 1;
    int strides[MAX_DIM] = {1};
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
        if (d > 0) {
            strides[d] = strides[d-1] * model->size[d-1];
        }
    }
    // Random site selection
    int site = rand() % total_size;
    int current_state = model->lattice[site];
    int new_state = rand() % model->q;
    // Calculate energy difference
    double delta_energy = 0.0;
    // For each dimension
    for (int d = 0; d < model->dim; d++) {
        // Positive neighbor
        int neighbor_pos = site + strides[d];
        // Handle boundary conditions
        if ((site / strides[d]) % model->size[d] == model->size[d] - 1) {
            if (model->boundary == 0) { // Periodic
                neighbor_pos -= model->size[d] * strides[d];
            } else { // Open
                neighbor_pos = -1;
            }
        }
        if (neighbor_pos >= 0 && neighbor_pos < total_size) {
            delta_energy -= (new_state == model->lattice[neighbor_pos]) ? -1.0 : 0.0;
            delta_energy -= (current_state == model->lattice[neighbor_pos]) ? 1.0 : 0.0;
        }
        // Negative neighbor
        int neighbor_neg = site - strides[d];
        // Handle boundary conditions
        if ((site / strides[d]) % model->size[d] == 0) {
            if (model->boundary == 0) { // Periodic
                neighbor_neg += model->size[d] * strides[d];
            } else { // Open
                neighbor_neg = -1;
            }
        }
        if (neighbor_neg >= 0 && neighbor_neg < total_size) {
            delta_energy -= (new_state == model->lattice[neighbor_neg]) ? -1.0 : 0.0;
            delta_energy -= (current_state == model->lattice[neighbor_neg]) ? 1.0 : 0.0;
        }
    }
    /// Metropolis acceptance criterion
    if (delta_energy <= 0.0 || (model->temperature > 0 && exp(-delta_energy / model->temperature) > (double)rand() / RAND_MAX)) {
        model->lattice[site] = new_state;
        model->energy += delta_energy;
        model->changed_sites[model->changed_count] = site;
        model->change_times[model->changed_count] = g_get_monotonic_time();
        model->changed_count++;
    }
    model->steps++;
    potts_calculate_magnetization(model);
}

/// Heat Bath algorithm step
void potts_heat_bath_step(PottsModel *model) {
    int total_size = 1;
    int strides[MAX_DIM] = {1};
    
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
        if (d > 0) {
            strides[d] = strides[d-1] * model->size[d-1];
        }
    }
    // Random site selection
    int site = rand() % total_size;
    int current_state = model->lattice[site];
    // Calculate probabilities for each possible state
    double probabilities[MAX_Q];
    double total = 0.0;
    for (int new_state = 0; new_state < model->q; new_state++) {
        double energy = 0.0;
        // For each dimension
        for (int d = 0; d < model->dim; d++) {
            // Positive neighbor
            int neighbor_pos = site + strides[d];
            // Handle boundary conditions
            if ((site / strides[d]) % model->size[d] == model->size[d] - 1) {
                if (model->boundary == 0) { // Periodic
                    neighbor_pos -= model->size[d] * strides[d];
                } else { // Open
                    neighbor_pos = -1;
                }
            }
            if (neighbor_pos >= 0 && neighbor_pos < total_size) {
                energy += (new_state == model->lattice[neighbor_pos]) ? -1.0 : 0.0;
            }
            // Negative neighbor
            int neighbor_neg = site - strides[d];
            // Handle boundary conditions
            if ((site / strides[d]) % model->size[d] == 0) {
                if (model->boundary == 0) { // Periodic
                    neighbor_neg += model->size[d] * strides[d];
                } else { // Open
                    neighbor_neg = -1;
                }
            }
            if (neighbor_neg >= 0 && neighbor_neg < total_size) {
                energy += (new_state == model->lattice[neighbor_neg]) ? -1.0 : 0.0;
            }
        }
        probabilities[new_state] = exp(-energy / model->temperature);
        total += probabilities[new_state];
    }
    // Normalize probabilities
    for (int i = 0; i < model->q; i++) {
        probabilities[i] /= total;
    }
    // Choose new state according to probabilities
    double r = (double)rand() / RAND_MAX;
    double cumulative = 0.0;
    int new_state = 0;
    for (; new_state < model->q - 1; new_state++) {
        cumulative += probabilities[new_state];
        if (r <= cumulative) break;
    }
    // Update lattice and recalculate energy if state changed
    if (new_state != current_state) {
        model->lattice[site] = new_state;
        potts_calculate_energy(model);
        model->changed_sites[model->changed_count] = site;
        model->change_times[model->changed_count] = g_get_monotonic_time();
        model->changed_count++;
    }
    model->steps++;
    potts_calculate_magnetization(model);
}

/// Swendsen-Wang cluster algorithm 
void potts_swendsen_wang_step(PottsModel *model, ClusterInfo *clusters) {
    int total_size = 1;
    int strides[MAX_DIM] = {1};
    for (int d = 0; d < model->dim; d++) {
        total_size *= model->size[d];
        if (d > 0) {
            strides[d] = strides[d-1] * model->size[d-1];
        }
    }
    // Reset cluster structures
    for (int i = 0; i < total_size; i++) {
        clusters->parent[i] = i;
        clusters->size[i] = 1;
    }
    // Bond activation step
    double p = 1.0 - exp(-1.0 / model->temperature);
    for (int i = 0; i < total_size; i++) {
        int current_state = model->lattice[i];
        // For each dimension (only positive neighbors to avoid double counting)
        for (int d = 0; d < model->dim; d++) {
            int neighbor = i + strides[d];
            // Handle boundary conditions
            if ((i / strides[d]) % model->size[d] == model->size[d] - 1) {
                if (model->boundary == 0) { // Periodic
                    neighbor -= model->size[d] * strides[d];
                } else { // Open
                    continue;
                }
            }
            if (neighbor < total_size && model->lattice[neighbor] == current_state) {
                if ((double)rand() / RAND_MAX < p) {
                    union_clusters(clusters->parent, clusters->size, i, neighbor);
                }
            }
        }
    }
    /// Label clusters
    clusters->total_clusters = 0;
    for (int i = 0; i < total_size; i++) {
        if (find_root(clusters->parent, i) == i) {
            clusters->labels[i] = clusters->total_clusters++;
        }
    }
    /// Assign cluster labels to all sites
    for (int i = 0; i < total_size; i++) {
        clusters->labels[i] = clusters->labels[find_root(clusters->parent, i)];
    }
    /// Flip clusters
    int *cluster_states = malloc(clusters->total_clusters * sizeof(int));
    for (int i = 0; i < clusters->total_clusters; i++) {
        cluster_states[i] = rand() % model->q;
    }
    /// Track changed sites
    model->changed_count = 0;
    for (int i = 0; i < total_size; i++) {
        int new_state = cluster_states[clusters->labels[i]];
        if (model->lattice[i] != new_state) {
            model->changed_sites[model->changed_count] = i;
            model->change_times[model->changed_count] = g_get_monotonic_time();
            model->changed_count++;
        }
        model->lattice[i] = new_state;
    }
    free(cluster_states);
    /// Recalculate energy and magnetization
    potts_calculate_energy(model);
    potts_calculate_magnetization(model);
    model->steps++;
}

/// Perform one simulation step
void potts_step(PottsModel *model, ClusterInfo *clusters) {
    model->changed_count = 0; // Reset changed sites
    switch (model->algorithm) {
        case 0: // Metropolis
            potts_metropolis_step(model);
            break;
        case 1: // Heat Bath
            potts_heat_bath_step(model);
            break;
        case 2: // Swendsen-Wang
            potts_swendsen_wang_step(model, clusters);
            break;
        default:
            potts_metropolis_step(model);
    }
}

void draw_callback(GtkDrawingArea *area, cairo_t *cr, int width, int height, gpointer user_data) {
    AppData *app = (AppData *)user_data;
    PottsModel *model = &app->model;
    if (!model->initialized) {
        cairo_set_source_rgb(cr, 1, 1, 1);
        cairo_paint(cr);
        return;
    }
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);
    gint64 current_time = g_get_monotonic_time();
    int rows = model->size[0];
    int cols = model->dim == 2 ? model->size[1] : model->size[2];
    int slice = model->dim == 3 ? model->size[0] / 2 : 0;
    double cell_width = (double)width / cols;
    double cell_height = (double)height / rows;
    double cell_size = MIN(cell_width, cell_height);
    GHashTable *changed_set = g_hash_table_new(NULL, NULL);
    for (int c = 0; c < model->changed_count; c++) {
        if ((current_time - model->change_times[c]) / 1000 < HIGHLIGHT_DURATION_MS) {
            g_hash_table_add(changed_set, GINT_TO_POINTER(model->changed_sites[c]));
        }
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int index = model->dim == 2 ? (i * cols + j) : (slice * rows * cols + i * cols + j);
            int state = model->lattice[index];
            GdkRGBA color = colors[state % MAX_Q];
            if (g_hash_table_contains(changed_set, GINT_TO_POINTER(index))) {
                cairo_set_source_rgba(cr, 
                    MIN(1.0, color.red + 0.3), 
                    MIN(1.0, color.green + 0.3), 
                    MIN(1.0, color.blue + 0.3), 
                    color.alpha);
            } else {
                cairo_set_source_rgba(cr, color.red, color.green, color.blue, color.alpha);
            }
            
            cairo_rectangle(cr, j * cell_size, i * cell_size, cell_size, cell_size);
            cairo_fill(cr);
            
            // Draw grid lines
            cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
            cairo_rectangle(cr, j * cell_size, i * cell_size, cell_size, cell_size);
            cairo_set_line_width(cr, 0.5);
            cairo_stroke(cr);
        }
    }
    
    g_hash_table_destroy(changed_set);
}

/// Simulation update function
gboolean update_simulation(gpointer user_data) {
    AppData *data = (AppData *)user_data;
    if (!data->running || data->paused || !data->model.initialized) {
        return G_SOURCE_CONTINUE;
    }
    potts_step(&data->model, &data->clusters);
    char energy_str[64];
    sprintf(energy_str, "Energy: %.2f", data->model.energy);
    gtk_label_set_text(GTK_LABEL(data->label_energy), energy_str);
    char mag_str[64];
    sprintf(mag_str, "Magnetization: %.4f", data->model.magnetization);
    gtk_label_set_text(GTK_LABEL(data->label_magnetization), mag_str);
    char steps_str[64];
    sprintf(steps_str, "Steps: %d", data->model.steps);
    gtk_label_set_text(GTK_LABEL(data->label_steps), steps_str);
    gtk_widget_queue_draw(data->drawing_area);
    return G_SOURCE_CONTINUE;
}

void start_simulation(GtkButton *button, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    if (!data->model.initialized) {
        data->model.q = gtk_range_get_value(GTK_RANGE(data->q_scale));
        data->model.dim = 2; // Fixed to 2D for mobile UI
        int size = gtk_range_get_value(GTK_RANGE(data->size_scale));
        for (int d = 0; d < data->model.dim; d++) {
            data->model.size[d] = size;
        }
        data->model.geometry = 0; // Square
        data->model.boundary = 0; // Periodic
        data->model.temperature = gtk_range_get_value(GTK_RANGE(data->temp_scale));
        data->model.algorithm = gtk_drop_down_get_selected(GTK_DROP_DOWN(data->combo_algorithm));
        potts_init(&data->model, &data->clusters);
    }
    if (!data->running) {
        data->running = true;
        data->paused = false;
        data->timeout_id = g_timeout_add(1000 / data->speed, update_simulation, data);
        gtk_widget_set_sensitive(data->start_button, FALSE);
        gtk_widget_set_sensitive(data->pause_button, TRUE);
        gtk_widget_set_sensitive(data->stop_button, TRUE);
        gtk_widget_set_sensitive(data->restart_button, TRUE);
        char energy_str[64];
        sprintf(energy_str, "Energy: %.2f", data->model.energy);
        gtk_label_set_text(GTK_LABEL(data->label_energy), energy_str);
        char mag_str[64];
        sprintf(mag_str, "Magnetization: %.4f", data->model.magnetization);
        gtk_label_set_text(GTK_LABEL(data->label_magnetization), mag_str);
        char steps_str[64];
        sprintf(steps_str, "Steps: %d", data->model.steps);
        gtk_label_set_text(GTK_LABEL(data->label_steps), steps_str);
        gtk_widget_queue_draw(data->drawing_area);
    }
}

void stop_simulation(GtkButton *button, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    if (data->timeout_id) {
        g_source_remove(data->timeout_id);
        data->timeout_id = 0;
    }
    data->running = false;
    data->paused = false;
    gtk_widget_set_sensitive(data->start_button, TRUE);
    gtk_widget_set_sensitive(data->pause_button, FALSE);
    gtk_widget_set_sensitive(data->stop_button, FALSE);
    gtk_button_set_label(GTK_BUTTON(data->pause_button), "Pause");
}

void pause_simulation(GtkButton *button, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    data->paused = !data->paused;
    gtk_button_set_label(button, data->paused ? "Resume" : "Pause");
}

void restart_simulation(GtkButton *button, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    if (data->timeout_id) {
        g_source_remove(data->timeout_id);
        data->timeout_id = 0;
    }
    data->model.q = gtk_range_get_value(GTK_RANGE(data->q_scale));
    data->model.dim = 2; // Fixed to 2D for mobile UI
    int size = gtk_range_get_value(GTK_RANGE(data->size_scale));
    for (int d = 0; d < data->model.dim; d++) {
        data->model.size[d] = size;
    }
    data->model.geometry = 0; // Square
    data->model.boundary = 0; // Periodic
    data->model.temperature = gtk_range_get_value(GTK_RANGE(data->temp_scale));
    data->model.algorithm = gtk_drop_down_get_selected(GTK_DROP_DOWN(data->combo_algorithm));
    potts_init(&data->model, &data->clusters);
    data->running = false;
    data->paused = false;
    gtk_widget_set_sensitive(data->start_button, TRUE);
    gtk_widget_set_sensitive(data->pause_button, FALSE);
    gtk_widget_set_sensitive(data->stop_button, FALSE);
    gtk_widget_set_sensitive(data->restart_button, TRUE);
    gtk_button_set_label(GTK_BUTTON(data->pause_button), "Pause");
    char energy_str[64];
    sprintf(energy_str, "Energy: %.2f", data->model.energy);
    gtk_label_set_text(GTK_LABEL(data->label_energy), energy_str);
    char mag_str[64];
    sprintf(mag_str, "Magnetization: %.4f", data->model.magnetization);
    gtk_label_set_text(GTK_LABEL(data->label_magnetization), mag_str);
    char steps_str[64];
    sprintf(steps_str, "Steps: %d", data->model.steps);
    gtk_label_set_text(GTK_LABEL(data->label_steps), steps_str);
    gtk_widget_queue_draw(data->drawing_area);
}

void speed_changed(GtkRange *range, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    data->speed = gtk_range_get_value(range);
    if (data->running && !data->paused) {
        if (data->timeout_id) {
            g_source_remove(data->timeout_id);
        }
        data->timeout_id = g_timeout_add(1000 / data->speed, update_simulation, data);
    }
}

void activate(GtkApplication *app, gpointer user_data) {
    AppData *data = (AppData *)user_data;
    data->window = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW(data->window), "Q-State Potts Model");
    gtk_window_set_default_size(GTK_WINDOW(data->window), 270, 480);
    gtk_window_set_resizable(GTK_WINDOW(data->window), FALSE);
    GtkCssProvider *provider = gtk_css_provider_new();
    gtk_css_provider_load_from_string(provider,
        "button {"
        "  padding: 10px;"
        "  border-radius: 8px;"
        "  font-size: 14px;"
        "  font-weight: bold;"
        "}"
        "scale {"
        "  margin: 8px 0;"
        "}"
        "label {"
        "  font-size: 12px;"
        "  margin-bottom: 4px;"
        "}"
        "combobox {"
        "  font-size: 12px;"
        "}");
    gtk_style_context_add_provider_for_display(gdk_display_get_default(),
        GTK_STYLE_PROVIDER(provider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
    g_object_unref(provider);
    GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
    gtk_window_set_child(GTK_WINDOW(data->window), box);
    gtk_widget_set_margin_start(box, 8);
    gtk_widget_set_margin_end(box, 8);
    gtk_widget_set_margin_top(box, 8);
    gtk_widget_set_margin_bottom(box, 8);
    data->drawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(data->drawing_area, 255, 336);
    gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(data->drawing_area), draw_callback, data, NULL);
    gtk_box_append(GTK_BOX(box), data->drawing_area);
    GtkWidget *control_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
    gtk_box_append(GTK_BOX(box), control_box);
    GtkWidget *q_label = gtk_label_new("Q (states):");
    gtk_box_append(GTK_BOX(control_box), q_label);
    data->q_scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 2, MAX_Q, 1);
    gtk_range_set_value(GTK_RANGE(data->q_scale), 2);
    gtk_box_append(GTK_BOX(control_box), data->q_scale);
    GtkWidget *size_label = gtk_label_new("Lattice Size:");
    gtk_box_append(GTK_BOX(control_box), size_label);
    data->size_scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 4, 100, 1);
    gtk_range_set_value(GTK_RANGE(data->size_scale), 50);
    gtk_box_append(GTK_BOX(control_box), data->size_scale);
    GtkWidget *temp_label = gtk_label_new("Temperature:");
    gtk_box_append(GTK_BOX(control_box), temp_label);
    data->temp_scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 0.01, 10.0, 0.01);
    gtk_range_set_value(GTK_RANGE(data->temp_scale), 1.0);
    gtk_box_append(GTK_BOX(control_box), data->temp_scale);
    GtkWidget *algo_label = gtk_label_new("Algorithm:");
    gtk_box_append(GTK_BOX(control_box), algo_label);
    data->combo_algorithm = gtk_drop_down_new_from_strings((const char *[]){"Metropolis", "Heat Bath", "Swendsen-Wang", NULL});
    gtk_box_append(GTK_BOX(control_box), data->combo_algorithm);
    GtkWidget *speed_label = gtk_label_new("Speed (FPS):");
    gtk_box_append(GTK_BOX(control_box), speed_label);
    data->speed_scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 1, 100, 1);
    gtk_range_set_value(GTK_RANGE(data->speed_scale), data->speed);
    gtk_box_append(GTK_BOX(control_box), data->speed_scale);
    g_signal_connect(data->speed_scale, "value-changed", G_CALLBACK(speed_changed), data);
    GtkWidget *button_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
    gtk_box_append(GTK_BOX(control_box), button_box);
    gtk_widget_set_halign(button_box, GTK_ALIGN_CENTER);
    data->start_button = gtk_button_new_with_label("Start");
    gtk_widget_set_size_request(data->start_button, 80, 48);
    gtk_box_append(GTK_BOX(button_box), data->start_button);
    g_signal_connect(data->start_button, "clicked", G_CALLBACK(start_simulation), data);
    data->pause_button = gtk_button_new_with_label("Pause");
    gtk_widget_set_size_request(data->pause_button, 80, 48);
    gtk_box_append(GTK_BOX(button_box), data->pause_button);
    gtk_widget_set_sensitive(data->pause_button, FALSE);
    g_signal_connect(data->pause_button, "clicked", G_CALLBACK(pause_simulation), data);
    data->stop_button = gtk_button_new_with_label("Stop");
    gtk_widget_set_size_request(data->stop_button, 80, 48);
    gtk_box_append(GTK_BOX(button_box), data->stop_button);
    gtk_widget_set_sensitive(data->stop_button, FALSE);
    g_signal_connect(data->stop_button, "clicked", G_CALLBACK(stop_simulation), data);
    data->restart_button = gtk_button_new_with_label("Restart");
    gtk_widget_set_size_request(data->restart_button, 80, 48);
    gtk_box_append(GTK_BOX(button_box), data->restart_button);
    g_signal_connect(data->restart_button, "clicked", G_CALLBACK(restart_simulation), data);
    data->label_energy = gtk_label_new("Energy: N/A");
    gtk_box_append(GTK_BOX(control_box), data->label_energy);
    data->label_magnetization = gtk_label_new("Magnetization: N/A");
    gtk_box_append(GTK_BOX(control_box), data->label_magnetization);
    data->label_steps = gtk_label_new("Steps: 0");
    gtk_box_append(GTK_BOX(control_box), data->label_steps);
    gtk_window_present(GTK_WINDOW(data->window));
}

int main(int argc, char **argv) {
    srand(time(NULL));
    init_colors();
    GtkApplication *app = gtk_application_new("org.example.potts", G_APPLICATION_DEFAULT_FLAGS);
    AppData data = {0};
    data.model.initialized = 0;
    data.clusters.labels = NULL;
    data.clusters.parent = NULL;
    data.clusters.size = NULL;
    data.running = false;
    data.paused = false;
    data.speed = 10.0;
    data.timeout_id = 0;
    g_signal_connect(app, "activate", G_CALLBACK(activate), &data);
    int status = g_application_run(G_APPLICATION(app), argc, argv);
    if (data.model.initialized) {
        free(data.model.lattice);
        free(data.model.changed_sites);
        free(data.model.change_times);
        free(data.clusters.labels);
        free(data.clusters.parent);
        free(data.clusters.size);
    }
    g_object_unref(app);
    return status;
}