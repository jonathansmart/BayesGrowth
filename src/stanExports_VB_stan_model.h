// Generated by rstantools.  Do not edit by hand.

/*
    BayesGrowth is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BayesGrowth is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BayesGrowth.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_VB_stan_model_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_VB_stan_model");
    reader.add_event(56, 54, "end", "model_VB_stan_model");
    return reader;
}
#include <stan_meta_header.hpp>
class model_VB_stan_model
  : public stan::model::model_base_crtp<model_VB_stan_model> {
private:
        int n;
        vector_d Age;
        vector_d Length;
        vector_d priors;
        vector_d priors_se;
public:
    model_VB_stan_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_VB_stan_model(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_VB_stan_model_namespace::model_VB_stan_model";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 1);
            current_statement_begin__ = 6;
            validate_non_negative_index("Age", "n", n);
            context__.validate_dims("data initialization", "Age", "vector_d", context__.to_vec(n));
            Age = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("Age");
            pos__ = 0;
            size_t Age_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < Age_j_1_max__; ++j_1__) {
                Age(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "Age", Age, 0);
            current_statement_begin__ = 7;
            validate_non_negative_index("Length", "n", n);
            context__.validate_dims("data initialization", "Length", "vector_d", context__.to_vec(n));
            Length = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("Length");
            pos__ = 0;
            size_t Length_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < Length_j_1_max__; ++j_1__) {
                Length(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "Length", Length, 0);
            current_statement_begin__ = 10;
            validate_non_negative_index("priors", "4", 4);
            context__.validate_dims("data initialization", "priors", "vector_d", context__.to_vec(4));
            priors = Eigen::Matrix<double, Eigen::Dynamic, 1>(4);
            vals_r__ = context__.vals_r("priors");
            pos__ = 0;
            size_t priors_j_1_max__ = 4;
            for (size_t j_1__ = 0; j_1__ < priors_j_1_max__; ++j_1__) {
                priors(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("priors_se", "2", 2);
            context__.validate_dims("data initialization", "priors_se", "vector_d", context__.to_vec(2));
            priors_se = Eigen::Matrix<double, Eigen::Dynamic, 1>(2);
            vals_r__ = context__.vals_r("priors_se");
            pos__ = 0;
            size_t priors_se_j_1_max__ = 2;
            for (size_t j_1__ = 0; j_1__ < priors_se_j_1_max__; ++j_1__) {
                priors_se(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "priors_se", priors_se, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 17;
            num_params_r__ += 1;
            current_statement_begin__ = 18;
            num_params_r__ += 1;
            current_statement_begin__ = 19;
            num_params_r__ += 1;
            current_statement_begin__ = 22;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_VB_stan_model() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 17;
        if (!(context__.contains_r("L0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable L0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("L0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "L0", "double", context__.to_vec());
        double L0(0);
        L0 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, L0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable L0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 18;
        if (!(context__.contains_r("Linf")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable Linf missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("Linf");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "Linf", "double", context__.to_vec());
        double Linf(0);
        Linf = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, Linf);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable Linf: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 19;
        if (!(context__.contains_r("k")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable k missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("k");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "k", "double", context__.to_vec());
        double k(0);
        k = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, k);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable k: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 22;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 17;
            local_scalar_t__ L0;
            (void) L0;  // dummy to suppress unused var warning
            if (jacobian__)
                L0 = in__.scalar_lb_constrain(0, lp__);
            else
                L0 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 18;
            local_scalar_t__ Linf;
            (void) Linf;  // dummy to suppress unused var warning
            if (jacobian__)
                Linf = in__.scalar_lb_constrain(0, lp__);
            else
                Linf = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 19;
            local_scalar_t__ k;
            (void) k;  // dummy to suppress unused var warning
            if (jacobian__)
                k = in__.scalar_lb_constrain(0, lp__);
            else
                k = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 22;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            // model body
            {
            current_statement_begin__ = 27;
            validate_non_negative_index("PredL", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> PredL(n);
            stan::math::initialize(PredL, DUMMY_VAR__);
            stan::math::fill(PredL, DUMMY_VAR__);
            current_statement_begin__ = 30;
            lp_accum__.add(normal_log<propto__>(Linf, get_base1(priors, 1, "priors", 1), get_base1(priors_se, 1, "priors_se", 1)));
            current_statement_begin__ = 31;
            lp_accum__.add(normal_log<propto__>(L0, get_base1(priors, 2, "priors", 1), get_base1(priors_se, 2, "priors_se", 1)));
            current_statement_begin__ = 32;
            lp_accum__.add(uniform_log<propto__>(k, 0, get_base1(priors, 3, "priors", 1)));
            current_statement_begin__ = 34;
            lp_accum__.add(uniform_log<propto__>(sigma, 0, get_base1(priors, 4, "priors", 1)));
            current_statement_begin__ = 38;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 39;
                stan::model::assign(PredL, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (Linf - ((Linf - L0) * stan::math::exp((-(k) * get_base1(Age, i, "Age", 1))))), 
                            "assigning variable PredL");
                current_statement_begin__ = 40;
                lp_accum__.add(normal_log(get_base1(Length, i, "Length", 1), get_base1(PredL, i, "PredL", 1), sigma));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("L0");
        names__.push_back("Linf");
        names__.push_back("k");
        names__.push_back("sigma");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_VB_stan_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double L0 = in__.scalar_lb_constrain(0);
        vars__.push_back(L0);
        double Linf = in__.scalar_lb_constrain(0);
        vars__.push_back(Linf);
        double k = in__.scalar_lb_constrain(0);
        vars__.push_back(k);
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 49;
            validate_non_negative_index("log_lik", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(n);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 50;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 51;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            normal_log(get_base1(Length, i, "Length", 1), (Linf - ((Linf - L0) * stan::math::exp((-(k) * get_base1(Age, i, "Age", 1))))), sigma), 
                            "assigning variable log_lik");
            }
            // validate, write generated quantities
            current_statement_begin__ = 49;
            size_t log_lik_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                vars__.push_back(log_lik(j_1__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_VB_stan_model";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "L0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Linf";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "k";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "L0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Linf";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "k";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_VB_stan_model_namespace::model_VB_stan_model stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
