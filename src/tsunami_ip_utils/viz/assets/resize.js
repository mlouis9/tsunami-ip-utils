window.addEventListener('resize', function() {
    const graphs = Array.from(document.querySelectorAll('.js-plotly-plot'));
    graphs.forEach(graph => {
        Plotly.Plots.resize(graph);
    });
});