#include "cgm.h"

// ind_from_zero - истина, если в файлах ig и jg содержат индексы с 0
cgm_solver::cgm_solver(const std::string& path, const bool& ind_from_zero)
{
    // n, max_iter и eps
    std::ifstream in(path + "kuslau.txt");

    in >> n >> max_iter >> eps;

    s.resize(n);
    az.resize(n);
    rr_vec.resize(n);
    z.resize(n);
    r.resize(n);

    in.close();

    // ig
    in.open(path + "ig.txt");
    ig.resize(n + 1);

    for (int i = 0; i < ig.size(); i++)
    {
        in >> ig[i];

        if (!ind_from_zero)
            --ig[i];
    }

    Ll.resize(ig.back());
    Ld.resize(n);

    in.close();

    // jg
    in.open(path + "jg.txt");
    jg.resize(ig.back());

    for (int i = 0; i < jg.size(); i++)
    {
        in >> jg[i];

        if (!ind_from_zero)
            --jg[i];
    }

    in.close();

    // ggl
    in.open(path + "ggl.txt");
    ggl.resize(ig.back());

    for (int i = 0; i < ggl.size(); i++)
        in >> ggl[i];

    in.close();

    // di
    in.open(path + "di.txt");
    diag.resize(n);

    for (int i = 0; i < diag.size(); i++)
        in >> diag[i];

    in.close();

    // pr
    in.open(path + "pr.txt");
    pr.resize(n);

    for (int i = 0; i < pr.size(); i++)
        in >> pr[i];

    pr_norm = norm(pr);
}

cgm_solver::cgm_solver(std::vector<int> &_ia, std::vector<int>& _ja,
    std::vector<double>& _di, std::vector<double>& _al, std::vector<double>& _b, double _eps)
{
    n = _di.size();
    max_iter = 1000;
    eps = _eps;
    //ggu = _al;
    ggl = _al;
    pr = _b;
    diag = _di;
    ig = _ia;
    jg = _ja;

    s.resize(n);
    az.resize(n);
    rr_vec.resize(n);
    z.resize(n);
    r.resize(n);
    Ll.resize(ig.back());
    Ld.resize(n);
    pr_norm = norm(pr);
}

int cgm_solver::get_n()
{
    return n;
}

int cgm_solver::get_total_iter()
{
    return total_iter;
}

// умножение матрицы СЛАУ на вектор x, результат помещается в вектор y
void cgm_solver::matrix_dot_vector(const std::vector<double>& x, std::vector<double>& y)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = diag[i] * x[i];

        for (int j = ig[i]; j < ig[i + 1]; j++)
        {
            y[i] += ggl[j] * x[jg[j]];
            y[jg[j]] += ggu[j] * x[i];
        }
    }
}

// разложение Холецкого для матрицы СЛАУ
void cgm_solver::llt()
{
    int p = 0, m = 0;

    for (int i = 0; i < n; i++)
    {
        Ld[i] = 0;

        for (int j = ig[i]; j < ig[i + 1]; j++)
        {
            Ll[j] = ggu[j];

            p = jg[j];
            m = ig[i];

            for (int k = ig[p]; k < ig[p + 1]; k++)
            {
                for (int q = m; q < j; q++)
                    if (jg[q] == jg[k])
                    {
                        Ll[j] -= Ll[q] * Ll[k];
                        m = q + 1;

                        break;
                    }
            }

            Ll[j] /= Ld[jg[j]];
            Ld[i] -= Ll[j] * Ll[j];
        }

        Ld[i] = sqrt(diag[i] + Ld[i]);
    }
}

// решение СЛАУ без предобуславливания
void cgm_solver::no_preconditioning(std::vector<double>& x0)
{
    // r = A * x0
    matrix_dot_vector(x0, r);

    // r = pr - A * x0 = pr - r
    vec_diff(pr, r, r);

    // z = r
    z = r;

    double r_norm = norm(r), rr = r_norm / pr_norm;

    while (rr >= eps && total_iter < max_iter)
    {
        // az = A * z
        matrix_dot_vector(z, az);

        // скалярное произведение p и r с предыдщуей итерации необходи-мо сохранить,
        // оно понадобится для вычисления bk после изменения p и r
        double scal_r_r = scalar_product(r, r), alpha_k = scal_r_r / scalar_product(az, z);

        for (int i = 0; i < n; i++)
        {
            x0[i] += alpha_k * z[i];
            r[i] -= alpha_k * az[i];
        }

        double beta_k = scalar_product(r, r) / scal_r_r;

        for (int i = 0; i < n; i++)
            z[i] = r[i] + beta_k * z[i];

        r_norm = norm(r);
        rr = r_norm / pr_norm;

        //std::cout << "iter #" << total_iter << " | rr = " << rr << std::endl;

        total_iter++;
    }

    std::cout << "RELATIVE RESIDUAL = " << relative_residual(x0) << std::endl;
    std::cout << "Iterations: " << total_iter << std::endl;

}

// решение СЛАУ с диагональным предобуславливанием
void cgm_solver::diag_preconditioning(std::vector<double>& x0)
{
    // r = A * x0
    matrix_dot_vector(x0, r);

    // r = pr - A * x0 = pr - r
    vec_diff(pr, r, r);

    // z = M^(-1) * r
    for (int i = 0; i < n; i++)
        z[i] = r[i] / diag[i];

    double r_norm = norm(r), rr = r_norm / pr_norm;

    while (rr >= eps && total_iter < max_iter)
    {
        // az = A * z
        matrix_dot_vector(z, az);

        double scal_m_inv_r = scalar_product(z, r), alpha_k = scal_m_inv_r / scalar_product(az, z);

        for (int i = 0; i < n; i++)
        {
            x0[i] += alpha_k * z[i];
            r[i] -= alpha_k * az[i];
        }

        // s = M^(-1) * r
        for (int i = 0; i < n; i++)
            s[i] = r[i] / diag[i];

        double beta_k = scalar_product(s, r) / scal_m_inv_r;

        for (int i = 0; i < n; i++)
            z[i] = s[i] + beta_k * z[i];

        r_norm = norm(r);
        rr = r_norm / pr_norm;

        //std::cout << "iter #" << total_iter << " | rr = " << rr << std::endl;

        total_iter++;
    }

    std::cout << "RELATIVE RESIDUAL = " << relative_residual(x0) << std::endl;
    std::cout << "Iterations: " << total_iter << std::endl;

}

// решение СЛАУ с предобуславливанием неполным разложением Холецкого
void cgm_solver::llt_preconditioning(std::vector<double>& x0)
{
    llt();
    if (x0.size() != r.size())
        x0.resize(r.size());
    // r = A * x0
    matrix_dot_vector(x0, r);

    // r = pr - A * x0 = pr - r
    vec_diff(pr, r, r);

    // z = M^(-1) * r
    solve_auxiliary_system(r, z);

    double r_norm = norm(r), rr = r_norm / pr_norm;

    while (rr >= eps && total_iter < max_iter)
    {
        // s = M^(-1) * r
        solve_auxiliary_system(r, s);

        // az = A * z
        matrix_dot_vector(z, az);

        double scal_m_inv_r = scalar_product(s, r), alpha_k = scal_m_inv_r / scalar_product(az, z);

        for (int i = 0; i < n; i++)
        {
            x0[i] += alpha_k * z[i];
            r[i] -= alpha_k * az[i];
        }

        // s = M^(-1) * r
        solve_auxiliary_system(r, s);

        double beta_k = scalar_product(s, r) / scal_m_inv_r;

        for (int i = 0; i < n; i++)
            z[i] = s[i] + beta_k * z[i];

        r_norm = norm(r);
        rr = r_norm / pr_norm;

        //std::cout << "iter #" << total_iter << " | rr = " << rr << std::endl;

        total_iter++;
    }

    std::cout << "RELATIVE RESIDUAL = " << relative_residual(x0) << std::endl;
    std::cout << "Iterations: " << total_iter << std::endl;

}

// решение вспомогательной СЛАУ (используется в методе с предобуславли-ванием разложением Холецкого)
void cgm_solver::solve_auxiliary_system(const std::vector<double>& f, std::vector<double>& x)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0.;

        for (int j = ig[i]; j < ig[i + 1]; j++)
            sum += Ll[j] * x[jg[j]];

        x[i] = (f[i] - sum) / Ld[i];
    }

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] /= Ld[i];

        for (int j = ig[i]; j < ig[i + 1]; j++)
            x[jg[j]] -= Ll[j] * x[i];
    }
}

double cgm_solver::relative_residual(const std::vector<double>& x)
{
    matrix_dot_vector(x, rr_vec);
    vec_diff(pr, rr_vec, rr_vec);

    return norm(rr_vec) / pr_norm;
}

void cgm_solver::vec_diff(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& res)
{
    for (int i = 0; i < x.size(); i++)
        res[i] = x[i] - y[i];
}

double cgm_solver::scalar_product(const std::vector<double>& x, const std::vector<double>& y)
{
    double sp = 0.;

    for (int i = 0; i < x.size(); i++)
        sp += x[i] * y[i];

    return sp;
}

double cgm_solver::norm(const std::vector<double>& x)
{
    return sqrt(scalar_product(x, x));
}

void cgm_solver::reset()
{
    total_iter = 0;

    s.assign(s.size(), 0.);
    az.assign(az.size(), 0.);
    rr_vec.assign(rr_vec.size(), 0.);
    z.assign(z.size(), 0.);
    r.assign(r.size(), 0.);

    Ld.assign(Ld.size(), 0.);
    Ll.assign(Ll.size(), 0.);
}

void cgm_solver::clear_all()
{
    total_iter = 0;

    ig.clear();
    ig.shrink_to_fit();
    jg.clear();
    jg.shrink_to_fit();
    ggu.clear();
    ggu.shrink_to_fit();
    diag.clear();
    diag.shrink_to_fit();
    pr.clear();
    pr.shrink_to_fit();
    r.clear();
    r.shrink_to_fit();
    z.clear();
    z.shrink_to_fit();
    rr_vec.clear();
    rr_vec.shrink_to_fit();
    az.clear();
    az.shrink_to_fit();
    s.clear();
    s.shrink_to_fit();
    Ll.clear();
    Ll.shrink_to_fit();
    Ld.clear();
    Ld.shrink_to_fit();
}