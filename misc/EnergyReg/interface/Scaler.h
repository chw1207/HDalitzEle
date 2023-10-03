#include <vector>
#include <algorithm>
#include <cmath>

#include "TObject.h"

class MinMaxScaler : public TObject {
private:
    std::vector<double> min_values_;
    std::vector<double> max_values_;

public:
    MinMaxScaler() {}

    void fit(const std::vector<std::vector<double>>& data) {
        // Initialize min and max values
        size_t num_features = data[0].size();
        min_values_.resize(num_features);
        max_values_.resize(num_features);
        std::fill(min_values_.begin(), min_values_.end(), std::numeric_limits<double>::max());
        std::fill(max_values_.begin(), max_values_.end(), -std::numeric_limits<double>::max());

        // Compute min and max values for each feature
        for (const auto& row : data) {
            for (size_t i = 0; i < num_features; i++) {
                min_values_[i] = std::min(min_values_[i], row[i]);
                max_values_[i] = std::max(max_values_[i], row[i]);
            }
        }
    }

    void transform(std::vector<double>& row) const {
        for (size_t i = 0; i < row.size(); i++) {
            double min_val = min_values_[i];
            double max_val = max_values_[i];
            row[i] = (row[i] - min_val) / (max_val - min_val);
        }
    }

    void inverse_transform(std::vector<double>& row) const {
        for (size_t i = 0; i < row.size(); i++) {
            double min_val = min_values_[i];
            double max_val = max_values_[i];
            row[i] = row[i] * (max_val - min_val) + min_val;
        }
    }

    ClassDef(MinMaxScaler, 1)
};


class StandardScaler : public TObject {
private:
    std::vector<double> mean_;
    std::vector<double> variance_;

public:
    StandardScaler() {}

    void fit(const std::vector<std::vector<double>>& data) {
        // Compute mean values for each feature
        size_t num_features = data[0].size();
        mean_.resize(num_features);
        std::fill(mean_.begin(), mean_.end(), 0.0);
        for (const auto& row : data) {
            for (size_t i = 0; i < num_features; i++) {
                mean_[i] += row[i];
            }
        }
        for (size_t i = 0; i < num_features; i++) {
            mean_[i] /= data.size();
        }

        // Compute variance values for each feature
        variance_.resize(num_features);
        std::fill(variance_.begin(), variance_.end(), 0.0);
        for (const auto& row : data) {
            for (size_t i = 0; i < num_features; i++) {
                variance_[i] += pow(row[i] - mean_[i], 2);
            }
        }
        for (size_t i = 0; i < num_features; i++) {
            variance_[i] /= data.size();
        }
    }

    void transform(std::vector<double>& row) const {
        for (size_t i = 0; i < row.size(); i++) {
            double mean_val = mean_[i];
            double variance_val = variance_[i];
            row[i] = (row[i] - mean_val) / sqrt(variance_val);
        }
    }

    void inverse_transform(std::vector<double>& row) const {
        for (size_t i = 0; i < row.size(); i++) {
            double mean_val = mean_[i];
            double variance_val = variance_[i];
            row[i] = row[i] * sqrt(variance_val) + mean_val;
        }
    }

    ClassDef(StandardScaler, 1)
};
