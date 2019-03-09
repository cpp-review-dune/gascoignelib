import sys
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib
matplotlib.use('pgf')

plt.rcParams.update({
    'text.usetex': True,
    'pgf.texsystem': 'pdflatex',
    'pgf.rcfonts': False,
    'font.family': 'serif',
    'font.serif': [],           # empty ones can be filled
    'font.sans-serif': [],
    'font.monospace': [],
    'font.size': 8,
    'figure.figsize': (3.15, 3),
    'pgf.preamble': [],
})


def bup_reader(file):
    if not file.is_file():
        return
    else:
        with file.open() as f:
            str = f.readline()
            dims = [int(s) for s in str.split() if s.isdigit()]
        bup_arr = np.genfromtxt(fname=file,
                                dtype=np.float64,
                                delimiter=' ',
                                skip_header=1,
                                usecols=range(dims[1]),
                                max_rows=dims[0])
        return bup_arr


def L2_error(path, approx_pattern, ref_pattern, out_path,
             pdf_export=True):
    app_list = list(sorted(path.glob(approx_pattern)))
    # exa_list = list(sorted(path.glob(exa_pattern)))
    ref_list = list(sorted(path.glob(ref_pattern)))
    # exa_length = len(exa_list)
    app_length = len(app_list)
    ref_length = len(ref_list)
    err_arr = np.zeros((ref_length, int(2 * app_length/ref_length)))
    # ref_err_arr = np.zeros((exa_length, 2))
    base = ['Pressure_L2_error_iteration_', 'Velocity_L2_error_iteration_']
    lbase = [r'$k=',
             r'$k=']
    cols = []
    [cols.append(b + str(i))
     for i in range(1, 1+int(app_length / ref_length)) for b in base]
    legend = []
    [legend.append(b + str(i) + '$')
     for i in range(1, 1+int(app_length / ref_length)) for b in lbase]
    cols_str = ' '.join(cols)
    
    for i in range(app_length):
        approx = bup_reader(app_list[i])
        ref = bup_reader(ref_list[i % ref_length])
        err = approx - ref
        v_err = np.linalg.norm(err[:, 1:2])/np.linalg.norm(ref[:, 1:2])
        p_err = np.linalg.norm(err[:, 0])/np.linalg.norm(ref[:, 0])
        err_arr[i % ref_length, 2 * int(i/ref_length)] = p_err
        err_arr[i % ref_length, 2 * int(i/ref_length) + 1] = v_err

    out_fname = pathlib.Path(out_path + '.txt')
    np.savetxt(out_fname, err_arr, delimiter=' ', header=cols_str, comments='')
    pdf_name = out_fname.with_suffix('.pdf')
    pgf_name = out_fname.with_suffix('.pgf')
    if pdf_export:
        with PdfPages(pdf_name.with_name('V' + pdf_name.name)) as pdf:
            X = np.linspace(1, ref_length, ref_length)
            plt.rc('text', usetex=True)
            plt.figure()
            plt.tight_layout(pad=0)
            plt.yscale('log')
            # plt.plot(X, ref_err_arr[:, 1], linewidth=.5, color='0.666')
            for i in range(1, 2 * int(app_length/ref_length), 2):
                plt.plot(X, err_arr[:, i], linewidth=.85, label=legend[i])
                plt.legend()
            pdf.savefig()
            plt.savefig(pgf_name.with_name('V' + pgf_name.name))
            plt.close()

        with PdfPages(pdf_name.with_name('P' + pdf_name.name)) as pdf:
            plt.rc('text', usetex=True)
            plt.figure()
            plt.tight_layout(pad=0)
            plt.yscale('log')
            # plt.plot(X, ref_err_arr[:, 0], linewidth=.5, color='0.666')
            for i in range(0, 2 * int(app_length/ref_length), 2):
                plt.plot(X, err_arr[:, i], linewidth=.85, label=legend[i])
                plt.legend()
            pdf.savefig()
            plt.savefig(pgf_name.with_name('P' + pgf_name.name))
            plt.close()
    else:
        pass


if __name__ == "__main__":
    curr_dir = pathlib.Path('.')
    if len(sys.argv) == 2:
        if sys.argv[1] == '2d':
            ref_pat = 'wm_2d/F001seq20/u.*.bup'
            approx_pat = 'wm_2d/Intervals20/C01tc51*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'wm_2d/C01tc51')

            approx_pat = 'wm_2d/Intervals20/C02tc52*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'wm_2d/C02tc52')

            approx_pat = 'wm_2d/Intervals20/C05tc55*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'wm_2d/C05tc55')

            approx_pat = 'wm_2d/Intervals20/C1tc6*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'wm_2d/C1tc6')
        elif sys.argv[1] == '3d':
            ref_pat = 'box3d/F001seq20/u.*.bup'
            approx_pat = 'box3d/Intervals20/C01tc51*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'box3d/C01tc51')

            approx_pat = 'box3d/Intervals20/C02tc52*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'box3d/C02tc52')

            approx_pat = 'box3d/Intervals20/C05tc55*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'box3d/C05tc55')

            approx_pat = 'box3d/Intervals20/C1tc6*/p.*.bup'
            L2_error(curr_dir, approx_pat, ref_pat, 'box3d/C1tc6')
    else:
        pass
