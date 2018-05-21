
   function test_make_f0_features()
   
   % need xmin, L, Nx, delv, vmax
   % load test file with parameters and expected f0 matrices
   test_list = {'spots','beam_vth0_s','beam_vth0_t','beam_vth0_x'};
   % add warm beams neutral, standing wave, and traveling wave
   
   for index = 1:numel(test_list)
       test = test_list{index};
       input_deck =  strcat( './unit_tests/make_f0vec_features_tests/', test);
       load(input_deck);
       f0_computed = make_f0_features(input_deck,spots,beams);
       
       % check that we are normalized
       delx = L/Nx;
       sum_f0_computed = delx*delv*sum(sum(f0_computed));
 
       
       if (max(max(f0 ~= f0_computed)))  || ~ ( (sum_f0_computed >.995) && (sum_f0_computed < 1.005) )
           strcat('Test f0 ', test, ' failed.')
       else
           strcat('Test f0 ', test, ' passes.')
       end
   end
  
   % save these into a parameter file to be loaded
   % call function with
   %     spot1-   - should get:
   %     spot2-   - should get:
   %     beam -  vth = 0, 's' -  should get:
   %     beam -  vth = 0, 't' -  should get:
   %     beam -  vth = 0, 'x' -  should get:
   %     beam -  vth =  , 's' -  should get:
   %     beam -  vth =  , 't' -  should get:
   %     beam -  vth =  , 'x' -  should get:
   end
   % 