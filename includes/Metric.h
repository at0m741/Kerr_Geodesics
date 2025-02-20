#pragma once	

#include "Geodesics.h"

class Metric {
	public:
		std::array<std::array<double, NDIM>, NDIM> gcov{};
		std::array<std::array<double, NDIM>, NDIM> gcon{};
		std::array<std::array<double, NDIM>, NDIM> gcovK{};
		std::array<std::array<double, NDIM>, NDIM> gconK{};
		std::array<std::array<double, NDIM>, NDIM> gconMinkowski{};
		std::array<std::array<double, NDIM>, NDIM> gcov_half{};
		std::array<std::array<double, NDIM>, NDIM> gcon_half{};
		std::array<std::array<double, NDIM>, NDIM> gcovMinkowski{};

		Metric() = default;

		~Metric() {}


		
		
		void calculate_metric(const std::array<double, NDIM>& x, 
				std::array<std::array<double, NDIM>, NDIM>& g,
				std::array<std::array<double, NDIM>, NDIM>& g_inv);
		void verify_metric(const std::array<std::array<double, NDIM>, NDIM>& g,
				const std::array<std::array<double, NDIM>, NDIM>& g_inv);
		void calculate_metric_kds(const std::array<double, NDIM>& x, 
				std::array<std::array<double, NDIM>, NDIM>& g,
				std::array<std::array<double, NDIM>, NDIM>& g_inv);
		void calculate_metric_kerr_newman(const std::array<double, NDIM>& x, 
				std::array<std::array<double, NDIM>, NDIM>& g,
				std::array<std::array<double, NDIM>, NDIM>& g_inv);

};
