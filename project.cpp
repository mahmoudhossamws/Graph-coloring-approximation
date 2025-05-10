#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>
#include <map>
#include <set>
#include <queue>
using namespace std;

struct GraphInstance {
    int id, n, m; // n = vertices, m = edges
    vector<vector<int>> adjacency_list;
    int optimal_colors, approx_colors;
    double ratio;
    string type;
};

unsigned int generateRandomSeed() {
    auto time_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    size_t addr_entropy = reinterpret_cast<size_t>(&time_seed);
    return static_cast<unsigned int>(time_seed ^ addr_entropy);
}

vector<GraphInstance> generateEdgeCases(int start_id = 1) {
    vector<GraphInstance> edge_cases;
    int id = start_id;
    
    // Edge Case 1: Complete Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Complete Graph";
        instance.n = 10;
        instance.m = 0; // Initialize edge count
        instance.m = instance.n * (instance.n - 1) / 2;
        instance.adjacency_list.resize(instance.n);
        
        for (int i = 0; i < instance.n; i++) {
            for (int j = 0; j < instance.n; j++) {
                if (i != j) {
                    instance.adjacency_list[i].push_back(j);
                }
            }
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 2: Bipartite Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Bipartite Graph";
        instance.n = 20;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        // Create a bipartite graph with two sets of 10 vertices each
        for (int i = 0; i < 10; i++) {
            for (int j = 10; j < 20; j++) {
                instance.adjacency_list[i].push_back(j);
                instance.adjacency_list[j].push_back(i);
                instance.m++;
            }
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 3: Star Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Star Graph";
        instance.n = 15;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        // Connect center (vertex 0) to all other vertices
        for (int i = 1; i < instance.n; i++) {
            instance.adjacency_list[0].push_back(i);
            instance.adjacency_list[i].push_back(0);
            instance.m++;
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 4: Cycle Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Cycle Graph";
        instance.n = 15;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        for (int i = 0; i < instance.n; i++) {
            instance.adjacency_list[i].push_back((i + 1) % instance.n);
            instance.adjacency_list[(i + 1) % instance.n].push_back(i);
            instance.m++;
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 5: Wheel Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Wheel Graph";
        instance.n = 15;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        // Connect center (vertex 0) to all other vertices
        for (int i = 1; i < instance.n; i++) {
            instance.adjacency_list[0].push_back(i);
            instance.adjacency_list[i].push_back(0);
            instance.m++;
        }
        
        // Connect the outer cycle
        for (int i = 1; i < instance.n; i++) {
            instance.adjacency_list[i].push_back(1 + (i % (instance.n - 1)));
            instance.m++;
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 6: Grid Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Grid Graph";
        int grid_size = 4; // 4x4 grid
        instance.n = grid_size * grid_size;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                int vertex = i * grid_size + j;
                
                // Connect to right neighbor
                if (j < grid_size - 1) {
                    instance.adjacency_list[vertex].push_back(vertex + 1);
                    instance.adjacency_list[vertex + 1].push_back(vertex);
                    instance.m++;
                }
                
                // Connect to bottom neighbor
                if (i < grid_size - 1) {
                    instance.adjacency_list[vertex].push_back(vertex + grid_size);
                    instance.adjacency_list[vertex + grid_size].push_back(vertex);
                    instance.m++;
                }
            }
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 7: Sparse Random Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Sparse Random Graph";
        instance.n = 17;  // Reduced from 25 to ensure n ≤ 25
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        mt19937 rng(generateRandomSeed());
        uniform_int_distribution<int> dist(0, instance.n - 1);
        
        // Add just enough edges to make it connected
        for (int i = 1; i < instance.n; i++) {
            int j = dist(rng) % i;
            instance.adjacency_list[i].push_back(j);
            instance.adjacency_list[j].push_back(i);
            instance.m++;
        }
        
        // Add a few more random edges
        for (int i = 0; i < 10; i++) {
            int u = dist(rng);
            int v = dist(rng);
            if (u != v && find(instance.adjacency_list[u].begin(), instance.adjacency_list[u].end(), v) == instance.adjacency_list[u].end()) {
                instance.adjacency_list[u].push_back(v);
                instance.adjacency_list[v].push_back(u);
                instance.m++;
            }
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 8: Dense Random Graph
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Dense Random Graph";
        instance.n = 17;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        mt19937 rng(generateRandomSeed());
        uniform_real_distribution<double> prob(0.0, 1.0);
        
        for (int i = 0; i < instance.n; i++) {
            for (int j = i + 1; j < instance.n; j++) {
                if (prob(rng) < 0.7) {  // 70% chance of an edge
                    instance.adjacency_list[i].push_back(j);
                    instance.adjacency_list[j].push_back(i);
                    instance.m++;
                }
            }
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 9: Clique Chain
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Clique Chain";
        int clique_size = 5;
        int num_cliques = 3;
        instance.n = clique_size * num_cliques;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        // Create separate cliques
        for (int c = 0; c < num_cliques; c++) {
            for (int i = 0; i < clique_size; i++) {
                for (int j = i + 1; j < clique_size; j++) {
                    int u = c * clique_size + i;
                    int v = c * clique_size + j;
                    instance.adjacency_list[u].push_back(v);
                    instance.adjacency_list[v].push_back(u);
                    instance.m++;
                }
            }
        }
        
        // Connect the cliques in a chain
        for (int c = 0; c < num_cliques - 1; c++) {
            int u = c * clique_size;
            int v = (c + 1) * clique_size;
            instance.adjacency_list[u].push_back(v);
            instance.adjacency_list[v].push_back(u);
            instance.m++;
        }
        edge_cases.push_back(instance);
    }
    
    // Edge Case 10: Odd Cycle Plus
    {
        GraphInstance instance;
        instance.id = id++;
        instance.type = "EDGE_CASE: Odd Cycle Plus";
        instance.n = 15;
        instance.m = 0; // Initialize edge count
        instance.adjacency_list.resize(instance.n);
        
        // Create a cycle of length 5
        for (int i = 0; i < 5; i++) {
            instance.adjacency_list[i].push_back((i + 1) % 5);
            instance.adjacency_list[(i + 1) % 5].push_back(i);
            instance.m++;
        }
        
        // Connect the rest in a star-like fashion
        for (int i = 5; i < instance.n; i++) {
            int j = i % 5;
            instance.adjacency_list[i].push_back(j);
            instance.adjacency_list[j].push_back(i);
            instance.m++;
        }
        edge_cases.push_back(instance);
    }
    
    return edge_cases;
}

vector<GraphInstance> generateInstances(int num_instances = 100) {
    mt19937 rng(generateRandomSeed());
    vector<GraphInstance> instances;
    
    struct InstanceType {
        int min_vertices, max_vertices;
        double edge_probability;
        string description;
    };
    
   vector<InstanceType> instance_types = {
    {5, 10, 0.3, "Small sparse graph"},
    {10, 17, 0.2, "Medium sparse graph"},
    {15, 17, 0.1, "Large sparse graph"},
    {5, 10, 0.7, "Small dense graph"},
    {10, 17, 0.5, "Medium dense graph"},
    {15, 17, 0.3, "Large dense graph"},
    {10, 17, 0.4, "Medium balanced graph"}
};
    
    for (int i = 0; i < num_instances; i++) {
        uniform_int_distribution<int> type_dist(0, instance_types.size() - 1);
        InstanceType type = instance_types[type_dist(rng)];
        
        GraphInstance instance;
        instance.id = i + 1;
        instance.type = "Random: " + type.description;
        instance.m = 0; // Initialize edge count
        
        uniform_int_distribution<int> n_dist(type.min_vertices, min(type.max_vertices, 17));
        instance.n = n_dist(rng);
        instance.adjacency_list.resize(instance.n);
        
        for (int u = 0; u < instance.n; u++) {
            for (int v = u + 1; v < instance.n; v++) {
                uniform_real_distribution<double> edge_dist(0.0, 1.0);
                if (edge_dist(rng) < type.edge_probability) {
                    instance.adjacency_list[u].push_back(v);
                    instance.adjacency_list[v].push_back(u);
                    instance.m++;
                }
            }
        }
        
        instances.push_back(instance);
    }
    
    return instances;
}

// Exact coloring using bitmask DP
int solveExactColoring(vector<vector<int>>& graph) {
    int n = graph.size();
    
    // For empty graph
    if (n == 0) return 0;
    
    // For very small graphs, use simple approach
    if (n <= 2) {
        if (n == 1) return 1;
        return (graph[0].size() > 0) ? 2 : 1;
    }
    
    // Check if graph is bipartite
    if (n <= 25) {
        vector<int> color(n, -1);
        bool is_bipartite = true;
        
        for (int start = 0; start < n && is_bipartite; start++) {
            if (color[start] != -1) continue;
            
            queue<int> q;
            q.push(start);
            color[start] = 0;
            
            while (!q.empty() && is_bipartite) {
                int u = q.front();
                q.pop();
                
                for (int v : graph[u]) {
                    if (color[v] == -1) {
                        color[v] = 1 - color[u];
                        q.push(v);
                    } else if (color[v] == color[u]) {
                        is_bipartite = false;
                        break;
                    }
                }
            }
        }
        
        if (is_bipartite) return 2;
    }
    
    // For larger graphs, use bitmask DP
    if (n > 25) {
        // Fallback to approximate for large graphs
        return n;
    }
    
    // Precompute adjacency bitmasks for faster checks
    vector<int> adj_mask(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j : graph[i]) {
            adj_mask[i] |= (1 << j);
        }
    }
    
    // DP[mask] = minimum number of colors needed to color the vertices in the mask
    vector<int> dp(1 << n, INT_MAX);
    dp[0] = 0; // Empty set needs 0 colors
    
    // Find all independent sets (potential color classes)
    vector<int> independent_sets;
    for (int mask = 1; mask < (1 << n); mask++) {
        bool is_independent = true;
        for (int i = 0; i < n && is_independent; i++) {
            if (mask & (1 << i)) {
                // Check if this vertex conflicts with any other vertex in the set
                if (mask & adj_mask[i]) {
                    is_independent = false;
                }
            }
        }
        if (is_independent) {
            independent_sets.push_back(mask);
        }
    }
    
    // DP recurrence
    for (int mask = 1; mask < (1 << n); mask++) {
        // Try to color the remaining uncolored vertices
        for (int indep : independent_sets) {
            // Check if this independent set is applicable
            if ((mask & indep) == indep) {
                int remaining = mask ^ indep;
                dp[mask] = min(dp[mask], dp[remaining] + 1);
            }
        }
    }
    
    return dp[(1 << n) - 1]; // All vertices colored
}

int solveApproximateColoring(vector<vector<int>>& graph) {
    int n = graph.size();
    
    // Special case: Empty graph
    if (n == 0) return 0;
    
    vector<int> color(n, -1);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int u, int v) {
        return graph[u].size() > graph[v].size();
    });

    for(int u : order){
        set<int> used;
        for (int v : graph[u])
            if (color[v] != -1) used.insert(color[v]);
        int c = 0;
        while (used.count(c)) c++;
        color[u] = c;
    }

    return *max_element(color.begin(), color.end()) + 1;
}

void saveToCSV(const vector<GraphInstance>& instances, const string& filename) {
    ofstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    
    file << "instance_id,type,num_vertices,num_edges,optimal_colors,approx_colors,ratio,max_degree" << endl;
    
    for (const auto& instance : instances) {
        int max_degree = 0;
        for (const auto& adj : instance.adjacency_list) {
            max_degree = max(max_degree, (int)adj.size());
        }
        
        file << instance.id << "," 
             << "\"" << instance.type << "\"" << "," 
             << instance.n << "," 
             << instance.m << "," 
             << instance.optimal_colors << "," 
             << instance.approx_colors << "," 
             << fixed << setprecision(4) << instance.ratio << ","
             << max_degree << endl;
    }
    
    file.close();
    cout << "Results saved to " << filename << endl;
}

void generateHistogram(const vector<GraphInstance>& instances) {
    map<double, int> histogram;
    double bin_width = 0.05;
    
    for (double bin = 1.0; bin <= 2.0; bin += bin_width) {
        histogram[bin] = 0;
    }
    
    for (const auto& instance : instances) {
        double bin = ceil(instance.ratio / bin_width) * bin_width;
        histogram[bin]++;
    }
    
    double max_ratio = 0.0;
    int worst_instance_id = -1;
    GraphInstance worst_instance;
    
    for (const auto& instance : instances) {
        if (instance.ratio > max_ratio) {
            max_ratio = instance.ratio;
            worst_instance_id = instance.id;
            worst_instance = instance;
        }
    }
    
    cout << "\n=== APPROXIMATION RATIO HISTOGRAM ===" << endl;
    cout << "Ratio Range   | Count | Visualization" << endl;
    cout << "----------------------------------------" << endl;
    
    for (const auto& [bin, count] : histogram) {
        if (count > 0) {
            cout << fixed << setprecision(2) 
                 << bin - bin_width << " - " << bin << " | "
                 << setw(5) << count << " | ";
            
            int bar_length = min(count, 50);
            for (int i = 0; i < bar_length; i++) cout << "#";
            
            if (count > 50) cout << "... (" << count << " total)";
            
            cout << endl;
        }
    }
    
    cout << "\n=== WORST CASE ANALYSIS ===" << endl;
    cout << "Maximum approximation ratio: " << fixed << setprecision(4) << max_ratio << endl;
    cout << "Found in Instance #" << worst_instance_id << " (" << worst_instance.type << ")" << endl;
    
    double avg_ratio = 0.0;
    for (const auto& instance : instances) avg_ratio += instance.ratio;
    avg_ratio /= instances.size();
    
    cout << "Average approximation ratio: " << fixed << setprecision(4) << avg_ratio << endl;
    
    map<string, vector<double>> type_ratios;
    for (const auto& instance : instances) {
        type_ratios[instance.type].push_back(instance.ratio);
    }
    
    cout << "\n=== PERFORMANCE BY INSTANCE TYPE ===" << endl;
    for (const auto& [type, ratios] : type_ratios) {
        double type_avg = 0.0;
        for (double r : ratios) type_avg += r;
        type_avg /= ratios.size();
        
        cout << "Type: " << type << "\n  Count: " << ratios.size() 
             << "\n  Avg Ratio: " << fixed << setprecision(4) << type_avg << endl;
    }
}

int main() {
    vector<GraphInstance> random_instances = generateInstances(10000);
    vector<GraphInstance> edge_instances = generateEdgeCases(random_instances.size() + 1);
    
    vector<GraphInstance> all_instances;
    all_instances.insert(all_instances.end(), random_instances.begin(), random_instances.end());
    all_instances.insert(all_instances.end(), edge_instances.begin(), edge_instances.end());
    
    cout << "Total instances: " << all_instances.size() << endl;
    int counter=0;
    for (auto& instance : all_instances) {
        cout<<counter<<endl;
        counter++;
        instance.optimal_colors = solveExactColoring(instance.adjacency_list);
        instance.approx_colors = solveApproximateColoring(instance.adjacency_list);
        instance.ratio = (instance.approx_colors > 0) ? 
                        static_cast<double>(instance.approx_colors) / instance.optimal_colors : 9999.9999;
    }
    
    saveToCSV(all_instances, "graph_coloring_results.csv");
    generateHistogram(all_instances);
    
    return 0;
}
