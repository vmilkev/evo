#include "model.hpp"
#include "solver_pcg.hpp"
#include "iointerface.hpp"
#include "parser.hpp"
#include "modelparser.hpp"
#include "effects.hpp"

#include <thread>
#include <map>

// test class definition
class Base {
public:
    void foo( std::map<size_t,double> &in_map, std::vector<size_t> &range_vect, size_t thr_id ) // non-static member function
    {
        std::cout << "Thread using non-static member function, with id "<<thr_id<< std::endl;
    }
    static void foo1() // static member function
    {
        std::cout << "Thread using static member function as callable" << std::endl;
    }
};

void foo( std::map<size_t,double> &in_map, std::map<size_t,double> &out_map, std::vector<size_t> &range_vect, size_t thr_id ) // non-static member function
{
    size_t i_first = 0;

    if ( thr_id != 0 )
        i_first = range_vect[thr_id - 1] + 1;

    size_t i_last = range_vect[thr_id];

    std::map<size_t,double>::iterator it = in_map.begin();
    std::advance(it, i_first);

    for(size_t i = 0; i <= (i_last - i_first); i++)
    {            
        out_map[ it->first ] = it->second;
        //std::cout<<"thread no. "<<thr_id<<", key & val: "<<it->first<<" "<<it->second<<"\n";
        it++;
    }
}

void foo2( std::map<size_t,double> &in_map, std::vector<std::map<size_t,double>> &out_vec, std::vector<size_t> &range_vect, size_t thr_id ) // non-static member function
{
    size_t i_first = 0;

    if ( thr_id != 0 )
        i_first = range_vect[thr_id - 1] + 1;

    size_t i_last = range_vect[thr_id];

    std::map<size_t,double>::iterator it = in_map.begin();
    std::advance(it, i_first);

    for(size_t i = 0; i <= (i_last - i_first); i++)
    {            
        out_vec[thr_id][ it->first ] = it->second;
        //std::cout<<"thread no. "<<thr_id<<", key & val: "<<it->first<<" "<<it->second<<"\n";
        it++;
    }
}

int main(void)
{
    try
    {
        // Testing std::map paralelization

        std::cout<<"Testing std::map paralelization ..."<<"\n";

        std::map<size_t, double> datamap;
        std::map<size_t, double> testmap;
        
        size_t n_records = 100;
        
        const auto processor_count = std::thread::hardware_concurrency(); //may return 0 when not able to detect

        size_t n_threads = processor_count;

        std::cout<<"processor count: "<<processor_count<<"\n";

        for(size_t i = 0; i < n_records; i++) // test data
        {
            if( i%2 == 0 )
                datamap[i+0] = (double)(i+0);
        }

        size_t map_size = datamap.size();

        size_t work_load = 0;
        work_load = (size_t)( map_size / n_threads );

        size_t max_load = 1; // max load per thread

        if (work_load <= max_load)
        {
            n_threads = (size_t)n_threads / 2.0;
            work_load = (size_t)( map_size / n_threads );

            if (work_load <= max_load)
            {
                n_threads = 1;
                work_load = map_size;
            }
        }

        std::vector<size_t> last_element;
        
        size_t last_index = work_load;
        for(size_t i = 0; i < n_threads - 1; i++)
        {
            std::cout<<"last index: "<<last_index<<"\n";

            last_element.push_back(last_index);
            last_index += work_load;
        }
        last_element.push_back(map_size-1);
        std::cout<<"last index: "<<map_size-1<<", n_threads: "<<n_threads<<"\n";
        
        /*std::vector<std::vector<size_t>> first_last( n_threads, std::vector<size_t>(2,0) );

        // need to correct last_element list so it include only last elements of a specific rows
        std::map<size_t,double>::iterator it = datamap.begin();
        size_t i_first, i_last;        
        //std::advance(it, i_first);
        //size_t i_row = it->first; // key of a very first element
        size_t count = 0;
        while(count < last_element.size())
        {
            size_t element = last_element[count];
            std::advance(it, element);
        }

        for(size_t i = 0; i < datamap.size(); i++)
        {
            size_t row = it->first;
            if ( row > last_element[count] )
            {
                if ()
                {
                    //
                }
                count++
            }
            it++;
        }*/

        /*Base b;
        std::thread th0(&Base::foo, &b, std::ref(datamap), std::ref(last_element), 0);
        std::thread th1(&Base::foo, &b, std::ref(datamap), std::ref(last_element), 1);
        std::thread th2(&Base::foo, &b, std::ref(datamap), std::ref(last_element), 2);
        std::thread th3(&Base::foo, &b, std::ref(datamap), std::ref(last_element), 3);
        std::thread th2(&Base::foo1);*/

        /*std::map<size_t, double> out0;
        std::map<size_t, double> out1;
        std::map<size_t, double> out2;
        std::map<size_t, double> out3;

        std::thread th0(foo, std::ref(datamap), std::ref(out0), std::ref(last_element), 0);
        std::thread th1(foo, std::ref(datamap), std::ref(out1), std::ref(last_element), 1);
        std::thread th2(foo, std::ref(datamap), std::ref(out2), std::ref(last_element), 2);
        std::thread th3(foo, std::ref(datamap), std::ref(out3), std::ref(last_element), 3);

        th0.join();
        th1.join();
        th2.join();
        th3.join();

        testmap.insert(out0.begin(), out0.end());
        testmap.insert(out1.begin(), out1.end());
        testmap.insert(out2.begin(), out2.end());
        testmap.insert(out3.begin(), out3.end());*/

        std::vector<std::thread> vec_threads;
        std::vector<std::map<size_t,double>> vec_maps;
        
        for(size_t i = 0; i < n_threads; i++)
        {
            std::map<size_t,double> out;
            vec_maps.emplace_back( out );
            //vec_maps.push_back( out );
        }
        for(size_t i = 0; i < n_threads; i++)
        {
            //vec_threads.emplace_back( foo2, std::ref(datamap), std::ref(vec_maps), std::ref(last_element), i );
            vec_threads.emplace_back( foo, std::ref(datamap), std::ref(vec_maps[i]), std::ref(last_element), i );
        }
        for (auto &v : vec_threads) {
            v.join();
            std::cout<<"next"<<std::endl;
        }
        for(size_t i = 0; i < n_threads; i++)
        {
            testmap.insert(vec_maps[i].begin(), vec_maps[i].end());
        }
        std::cout<<"finish"<<std::endl;

        for (auto const &[key, val]: datamap)
        {
            std::cout<<"DATA key, val: "<<key<<" "<<val<<"\n";
        }

        for (auto const &[key, val]: testmap)
        {
            std::cout<<"TEST key, val: "<<key<<" "<<val<<"\n";
        }

        std::cout<<"Completed."<<"\n";

        exit(1);
        // -------------------------
        std::cout << "Parser:"
                  << "\n";

        //evolm::Modelparser m;

        //m.eval_expr( "y1,y2 , another_var, ,= x1+ x2 - 1 +(1+f|x3*x4) - x4 +(-1|x4) + x1+x2:x3 + x4*x5*x6 - x4:x5:x6" );

        //exit(0);
        // --------------------------
        std::cout << "Selected file reading:"
                  << "\n";

        /*
            Example:
            id   obs_f var_i1 var_f2 var_cat var_str var_str2
            1002 12.2  20     51.1   1       aple    true
            1003 15.5  30     10     2       plum    true
            1004 21.0  45     562    3       aple    false
            1052 30.5  50     452    3       plum    true
            1062 40    61     231    4       tomato  false
            1072 51.3  71     125    2       tomato  true
            1082 -9.0  80     121    1       plum    true
            1092 70.01 91     121    1       aple    false
            1102 82.12 10     110.0  4       tomato  false
        */
        evolm::IOInterface in;
        evolm::Modelparser m;

        in.set_fname("tests/data/diverse_data.dat");
        in.set_missing(-9.0);

        m.eval_expr("obs_f = -1 + var_cat + var_i1 + var_str*var_str2 + id");
        
        std::vector<std::string> obs;
        std::vector<std::string> fixed;
        std::vector<std::string> random;

        m.get_modelvars(obs,fixed,random);

         evolm::Effects obs_res;
        
        in.fgetvar(obs[0], obs_res);
        obs_res.print( obs[0] );

        // here we still need to parse "." (each col elementwise multiplication)
        // and concatenate effect matrices (+);
        // and vector of intercept (1).
        for (size_t i = 0; i < fixed.size(); i++)
        {
             evolm::Effects fix_res;
            in.fgetvar(fixed[i], fix_res);
            fix_res.print( fixed[i] );
        }

        // -------------------------
        exit(0);

        // Testing model

        evolm::Pcg solver;
        evolm::Model model;

        std::vector<float> iR{0.041}; //, 0.041, 0.041, 0.041};

        std::vector<float> iG1{10.1}; //, 0.1, 0.1, 0.1};

        // size_t dim = 1;

        std::cout << "In main."
                  << "\n";

        model.append_residual(iR, 1);

        model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0

        size_t type1 = 0;
        float type2 = 0.0f;

        model.append_effect("tests/data/model_4/obs_489_snp_1000.txt", type1); // eff := 0
        model.append_effect("tests/data/model_4/fixed_1.dat", type2);          // eff := 1

        std::vector<int> corr_eff{0};

        std::vector<int> eff_trate{1, 0};
        int obs_trate = 0;

        std::string identity("I");

        model.append_corrstruct(iG1, 1, identity, 1000, corr_eff);

        model.append_traitstruct(obs_trate, eff_trate);
        // model.append_traitstruct(obs_trate, eff_trate);

        solver.append_model(model);

        solver.solve(2);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}