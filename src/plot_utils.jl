"""
    meshgrid(xs, ys)

Create a meshgrid from the x and y coordinates specified by `xs`` and `ys`.`
"""
function meshgrid(xs, ys)
    Xs = [x for x in xs for y in ys]
    Ys = [y for x in xs for y in ys]

    return Xs, Ys
end

"""
    mesh_vector_field(Xs, Ys, f::Function, scale::Float64, normalize_arrows::Bool)

Compute the vector field over the mesh specified by `Xs` and `Ys`.

If the norm of the vector field is zero anywhere then its likely we have some NaNs in dx and
dy. If so, just replace any `NaN` values with zero.
"""
function mesh_vector_field(Xs, Ys, f::Function, scale::Float64, normalize_arrows::Bool)
    # Check if we should normalize all vectors to unit length and then scale them up/down
    if normalize_arrows
        dfs = scale * normalize.(f.(Xs, Ys))
    else
        dfs = scale * f.(Xs, Ys)
    end

    # Reshape into individual vectors for the x and y components of vector field
    dx = [df[1] for df in dfs]
    dy = [df[2] for df in dfs]

    # Check if we divided by zero anywhere
    replace!(dx, NaN => 0.0)
    replace!(dy, NaN => 0.0)

    return dx, dy
end

"""
    vector_field_colors(Xs, Ys, f::Function)

Map the magnitude of a vector to a color.

I don't know why this works. The answer is taken from
https://discourse.julialang.org/t/quiver-plot-plots-pyplot-with-different-colors-depending-on-the-length-of-the-arrow/59577/5
"""
function vector_field_colors(Xs, Ys, f::Function)
    c = norm.(f.(Xs, Ys))
    c = [c c]'
    c = repeat([c...], inner=2)

    return c
end

# Make some cool colors
myblue = "0x0772B3"
myred = "0xF0615C"
mygreen = "0x009F73"
mypurple = "0x786EB3"
myyellow = "0xE7A122"
mycyan = "0x5DB4E5"
mypink = "0xD95BA1"
mypalette = [myblue, myred, mygreen, mypurple, myyellow, mycyan, mypink]
mycolors = [
    parse.(Colorant, mypalette)[1],
    Colors.JULIA_LOGO_COLORS[1],
    Colors.JULIA_LOGO_COLORS[2],
    Colors.JULIA_LOGO_COLORS[4],
    parse.(Colorant, mypalette)[5],
]
get_colors() = mycolors
get_color_palette() =  parse.(Colorant, mypalette)

# Default plot settings for PGFplots
ratio = 1.15
fig_width = 4.0
fig_height = fig_width/1.15
ax_theme = @pgf {
    "semithick", 
    tick_style={"semithick", color="black"}, 
    legend_style={draw="none", fill="none", legend_cell_align="left", font="\\footnotesize"},
    max_space_between_ticks=100,
    width="$fig_width"*"in",
    height="$fig_width"*"in",
    tick_label_style={font="\\footnotesize"},
    xlabel_style={font="\\small"},
    ylabel_style={font="\\small"},
    xlabel_shift="-4pt",
    ylabel_shift="-4pt",
}
plt_theme = @pgf {"very thick", "no marks", line_join="round", line_cap="round"}
get_ax_theme() = ax_theme
get_plt_theme() = plt_theme
plot_defaults() = get_ax_theme(), get_plt_theme(), get_colors()

# Plot settings for IEEE style figures
ieee_textwidth = 7.0
ieee_column_sep = 0.2
# ieee_column_width = (ieee_textwidth - ieee_column_sep)/2
ieee_column_width = ieee_textwidth/2
ieee_ratio = 1.5 # TODO: determine if we want something else, maybe 1.5

ieee_theme = @pgf {
    "semithick", 
    tick_style={"semithick", color="black"}, 
    legend_style={draw="none", fill="none", legend_cell_align="left", font="\\footnotesize"},
    max_space_between_ticks=100,
    tick_label_style={font="\\footnotesize"},
    xlabel_style={font="\\small"},
    ylabel_style={font="\\small"},
    xlabel_shift="-4pt",
    ylabel_shift="-4pt",
}

# Setting for large figures - take up an entire column
ieee_width_large = ieee_column_width
ieee_height_large = ieee_width_large/ieee_ratio
ieee_theme_large = @pgf {
    ieee_theme...,
    width="$ieee_width_large"*"in",
    height="$ieee_height_large"*"in",
}

# Setting for medium figures - take up half column
ieee_width_medium = ieee_textwidth/3.0
ieee_height_medium = ieee_width_medium/ieee_ratio
ieee_theme_medium = @pgf {
    ieee_theme...,
    width="$ieee_width_medium"*"in",
    height="$ieee_height_medium"*"in",
}

# Setting for small figures - used for two side by side figures in one column
ieee_width_small = ieee_column_width/2
ieee_height_small = ieee_width_small/ieee_ratio
ieee_theme_small = @pgf {
    ieee_theme...,
    width="$ieee_width_small"*"in",
    height="$ieee_height_small"*"in",
}

# Function to get IEEE themes
get_ieee_theme_large() = ieee_theme_large
get_ieee_theme_medium() = ieee_theme_medium
get_ieee_theme_small() = ieee_theme_small

# Place Text in plots
TextNode(x, y, text) = [raw"\node at ", Coordinate(x, y), text, ";"]