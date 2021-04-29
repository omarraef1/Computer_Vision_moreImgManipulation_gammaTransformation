
function write_figure(base_fname)
    figure_dir = 'figures/';
    have_epstopdf = 1;

    if (strlength(figure_dir) > 0) 
        cmd = sprintf('mkdir -p %s', figure_dir');
        system(cmd);
    end 

    eps_fname = sprintf('%s%s.eps', figure_dir, base_fname);
    print('-depsc2', eps_fname)

    if (have_epstopdf)
        pdf_fname = sprintf('%s%s.pdf', figure_dir, base_fname);
        cmd = sprintf('epstopdf %s --outfile="%s"', eps_fname, pdf_fname);
        system(cmd);
    end
end

