#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <iostream>
#include <random>
#include <limits>
using namespace std;

struct Edge {
    int to;
    vector<double> objectives;

    Edge() : to(0), objectives({}) {}
    Edge(int to_, vector<double> obj) : to(to_), objectives(obj) {}
};

vector<vector<Edge>> g;
int n, m, d, source, destination;
const int MAX_ITER = 1000;

double temperature = 1000.0;
double cooling_rate = 0.995;

struct Solution {
    vector<int> path;
    vector<double> objectives;

    Solution() {}
    Solution(vector<int> p, vector<double> obj) : path(p), objectives(obj) {}
};

bool dominates(const vector<double>& a, const vector<double>& b) {
    bool better_in_any = false;
    for (int i = 0; i < d; ++i) {
        if (a[i] > b[i]) return false;
        if (a[i] < b[i]) better_in_any = true;
    }
    return better_in_any;
}

vector<double> calculateObjectives(const vector<int>& path) {
    vector<double> objectives(d, 0.0);
    for (size_t i = 1; i < path.size(); ++i) {
        int from = path[i - 1];
        int to = path[i];
        for (const Edge& edge : g[from]) {
            if (edge.to == to) {
                for (int j = 0; j < d; ++j) {
                    objectives[j] += edge.objectives[j];
                }
                break;
            }
        }
    }
    return objectives;
}

double sumObjectives(const vector<double>& objectives) {
    double sum = 0.0;
    for (double obj : objectives) {
        sum += obj;
    }
    return sum;
}

vector<int> generateRandomPath() {
    vector<int> path;
    vector<bool> visited(n, false);
    int current = source;
    path.push_back(current);
    visited[current] = true;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, n - 1);

    while (current != destination) {
        vector<int> neighbors;
        for (const Edge& edge : g[current]) {
            if (!visited[edge.to]) {
                neighbors.push_back(edge.to);
            }
        }
        if (neighbors.empty()) break; 
        current = neighbors[dis(gen) % neighbors.size()];
        path.push_back(current);
        visited[current] = true;
    }
    return path;
}

Solution simulatedAnnealing(Solution initial_solution) {
    Solution current_solution = initial_solution;
    Solution best_solution = initial_solution;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    while (temperature > 1.0) {
        Solution new_solution = current_solution;

        int swap_idx1 = dis(gen) * (new_solution.path.size() - 1) + 1;
        int swap_idx2 = dis(gen) * (new_solution.path.size() - 1) + 1;
        swap(new_solution.path[swap_idx1], new_solution.path[swap_idx2]);

        new_solution.objectives = calculateObjectives(new_solution.path);

        double current_cost = sumObjectives(current_solution.objectives);
        double new_cost = sumObjectives(new_solution.objectives);

        if (new_cost < current_cost || exp((current_cost - new_cost) / temperature) > dis(gen)) {
            current_solution = new_solution;
        }

        if (sumObjectives(current_solution.objectives) < sumObjectives(best_solution.objectives)) {
            best_solution = current_solution;
        }

        temperature *= cooling_rate;
    }

    return best_solution;
}

Solution antLionOptimizer() {
    vector<Solution> population;
    for (int i = 0; i < n; ++i) {
        vector<int> path = generateRandomPath();
        population.emplace_back(path, calculateObjectives(path));
    }

    Solution best_solution = population[0];
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        for (Solution& sol : population) {
            sol = simulatedAnnealing(sol);
            if (sumObjectives(sol.objectives) < sumObjectives(best_solution.objectives)) {
                best_solution = sol;
            }
        }
    }
    return best_solution;
}

int main() {
    fstream cin("input.txt");
    ofstream cout("output.txt");

    cin >> n >> m >> d;
    cin >> source >> destination;
    g.resize(n);

    for (int i = 0; i < m; i++) {
        int from, to;
        cin >> from >> to;
        vector<double> temp(d);
        for (int j = 0; j < d; j++) {
            cin >> temp[j];
        }
        g[from].emplace_back(to, temp);
        g[to].emplace_back(from, temp); 
    }

    Solution optimal_solution = antLionOptimizer();

    if (!optimal_solution.path.empty()) {
        for (int v : optimal_solution.path) {
            cout << v << " ";
        }
        cout << endl;
    }
    else {
        cout << "No path found" << endl;
    }

    return 0;
}
