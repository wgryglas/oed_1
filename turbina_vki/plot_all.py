

if __name__ == "__main__":
    import settings
    import plot_comparison_on_boundary
    import plot_modes_over_boundary
    import plot_sepctrum
    import plot_optimized_points

    settings.run_modules(plot_optimized_points, plot_comparison_on_boundary, plot_modes_over_boundary, plot_sepctrum)
