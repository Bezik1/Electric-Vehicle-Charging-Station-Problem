#include "Maps.hpp"

Maps::Maps(
    int grid_width,
    int grid_height,
    const std::vector<std::vector<double>>& distances_costs_map,
    const std::vector<std::vector<int>>& poi_map,
    const std::vector<std::vector<double>>& demand_map,
    const std::vector<std::vector<double>>& land_rental_cost_map
) : grid_width{grid_width},
    grid_height{grid_height},
    distances_costs_map{distances_costs_map},
    poi_map{poi_map},
    demand_map{demand_map},
    land_rental_cost_map{land_rental_cost_map}
{
    distances_from_to.assign(grid_width, std::vector<std::vector<std::vector<double>>>(
        grid_height, 
        std::vector<std::vector<double>>(grid_width, std::vector<double>(
                                                        grid_height, std::numeric_limits<double>::infinity()))
    ));

    for (int m = 0; m < grid_width; ++m)
        for (int n = 0; n < grid_height; ++n)
        {
            if (demand_map[m][n] <= 0) continue;

            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                    distances_from_to[m][n][i][j] = calculate_distance(m, n, i, j);
            }
        }
}

std::pair<double, std::pair<int, int>> Maps::find_nearest_valid(int start_x, int start_y) const
{
    double distance_to_nearest_valid = 0;

    if (distances_costs_map[start_x][start_y] != INF_VAL)
        return { distance_to_nearest_valid, {start_x, start_y} };

    std::queue<std::pair<int, int>> q;
    std::vector<std::vector<bool>> visited(grid_width, std::vector<bool>(grid_height, false));

    q.push({start_x, start_y});
    visited[start_x][start_y] = true;

    while(!q.empty())
    {
        distance_to_nearest_valid += 1;
        auto [cx, cy] = q.front();
        q.pop();

        for(const auto& [dx, dy] : DIRECTIONS)
        {
            int nx = cx + dx;
            int ny = cy + dy;

            if(nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height && !visited[nx][ny])
            {
                if (distances_costs_map[nx][ny] != INF_VAL)
                    return { distance_to_nearest_valid, {nx, ny} };

                visited[nx][ny] = true;
                q.push({nx, ny});
            }
        }
    }
    return { distance_to_nearest_valid, {-1, -1} };
}

double Maps::calculate_distance(int m, int n, int i, int j) const
{
    using Node = std::pair<double, std::pair<int, int>>;

    if (m == i && n == j)
        return 0.0;

    if (distances_costs_map[i][j] == INF_VAL)
        return std::numeric_limits<double>::infinity();

    auto [extra_d, actual_start] = find_nearest_valid(m, n);
    if (actual_start.first == -1)
        return std::numeric_limits<double>::infinity();

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
    std::vector<std::vector<double>> dist(
        grid_width, 
        std::vector<double>(grid_height, std::numeric_limits<double>::infinity())
    );

    dist[actual_start.first][actual_start.second] = extra_d;
    queue.push({extra_d, actual_start});

    while (!queue.empty())
    {
        auto [current_d, coord] = queue.top();
        auto [x, y] = coord;
        queue.pop();

        if (current_d > dist[x][y]) continue;
        if (x == i && y == j) return current_d;

        for (const auto& [dx, dy] : DIRECTIONS)
        {
            int nx = x + dx; int ny = y + dy;

            if (nx < 0 || nx >= grid_width || ny < 0 || ny >= grid_height || !isPassable(nx, ny))
                continue;

            double weight = (distances_costs_map[nx][ny] <= 0) ? 1.0 : distances_costs_map[nx][ny];
            if (current_d + weight < dist[nx][ny]) {
                dist[nx][ny] = current_d + weight;
                queue.push({dist[nx][ny], {nx, ny}});
            }
        }
    }
    return std::numeric_limits<double>::infinity();
}