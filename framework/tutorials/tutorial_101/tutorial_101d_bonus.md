## Part d: bonus tasks

### (1) Error handling - Part (i): Using the `*.lst` file

There are two different kinds of errors a model can have: (i) errors in the structure of the model; and (ii) errors deriving from infeasibilities which make the model not solvable. In this first bonus task, we look at category (i).
Whenever there is talk about `*.dat` files in the following, they could also be `*.csv` files if you chose the option `"fileformat="csv"` when calling the Instance.write() method.

#### (1.1) Running a model with a mistake in a `*.dat` or `*.csv` file

It is handy to have an idea of how to find and handle error messages in REMix. We will produce such an error in this first bonus task. To do so, change the word "lifeTime" in line 234 to "lifetime", for example. Afterwards, run the model again. It won't solve but give you a `*** Status: Compilation error(s)` message in the log which is terminal output. Find that message.

#### (1.2) Finding the error

The status message is not very helpful in figuring out what might have gone wrong with your model. The `*.lst` file can be, however. Access the `*.lst` file as described above when introducing how to run the GAMS model. You can try out both mentioned possibilities.

#### (1.3) Understanding the error 1

If you caused the error proposed in task (1.1), you will find the following error in the `*.lst` file:

```
1714  ,,lifeTime,activityUpperLimit,lifetime
****                                       $172
**** LINE      1 INCLUDE     C:\work\remix\framework\225d\/converter_techparam.csv
**** LINE     79 INCLUDE     C:\work\remix\framework\/model\/core\/converter.gms
**** LINE     24 INCLUDE     C:\work\remix\framework\/model\/framework\/remix.gms
**** LINE     22 INPUT       C:\work\remix\framework\run_remix.gms
**** 172  Element is redefined
```

The position of the dollar sing ("$") followed by a number will always tell you where the error occurred. In this case, the solver does not seem to be happy about the word "lifetime" and tells you (in the bottom-most line) that an "Element is redefined" (which is the description of error number 172).

You can also see in which file the error occurred. Here, it is in `converter_techparam.csv` (during the GAMS call, the `*.dat` files we created are converted to `*.csv` files). This will in some cases already help you to understand your mistake.

#### (1.4) Understanding the error 2

If the method presented in task (1.3) does not suffice to understand the error in the `*.dat` file, this might be because you are not sure how this file is supposed to look like correctly. Luckily, you can easily find that out by checking a correctly set-up model.

You can find some in the folder `testing`, which includes most basic `*.dat` files for different functions for testing purposes. This folder is situated on the same layer as the `tutorials` folder.

If you now go to, e.g., `testing/instances/minimal_lp/data`, you will find a lot of `*.dat` files, one of which called `converter_techParam.dat`. Open it and compare it to the `converter_techParam.dat` in your data_dir. Now, you can see that there is a "lifetime" column too much in your `*.dat` file and realize that the "t" should be an upper-case letter.

### (2) Error handling - Part (ii): Using the log (infeasibilities)

In this second bonus task, we will have a look at infeasibilities and see how we can find and understand them.

#### (2.1) Causing an infeasibility and finding the error message

In order to see an example of a model that is infeasible to solve, set the upper limit for the CCGT plant to 1 GW.
After you have done that, run the model again and look for the infeasibility message in the terminal output. It should read:

`Row 'c16' infeasible, all entries at implied bounds.`

This is mighty unhelpful, you might think. We think that anyway. However, there is a way to make this message a bit more meaningful to human readers.

N.B.: In case you have wondered: you will also get a "KeyError" if you run the full script. This error happens in the evaluation part of the tutorial, however, and is just a consequence of the infeasibility caused above.

#### (2.2) Understanding the infeasibility

We can get a much more meaningful infeasibility message by changing just one thing, which is by including `names=1` as command line option to the GAMS function call. Include this argument in line 545, run the model again and look again for the infeasibility message. It should now read:

`Row 'Eq_converter_activityUpperLimit(t0001,DE_model,2030,CCGT,2030)' infeasible, all entries at implied bounds.`

Nice! This is something we can work with! So, the problem is in the equation for the activityUpperLimit of a converter. If we read on, we can read that this converter is the CCGT plant, which is no surprise, since we have just changed its maximum installable capacity to 1 GW, which does not seem to be enough.
We can also learn that the infeasibility already happens in the first time step, in the model region "DE_model" for the optimization year 2030.

#### (2.3) Giving out a dedicated log file

You might have wondered why the `*.lst` file that we used in bonus task 1 is not helpful at all here. It will only tell you about a "normal completion" of the model run, and apart from that only that the model status is "infeasible".
If you find the terminal output unsufficient to work with, there is the possibility to get a dedicated log file that will show the infeasibility message. This, again, can be done using a command line option.

If you would like REMix to store a `run_remix.log` file for you, just change the command line argument `lo=3` to `lo=4`. After another model run, you will find the log file in the same folder as the `*.lst` file. Its contents should already be familiar to you from the terminal outputs of the previous model runs.

### (3) Adding a new technology for wind offshore

In addition to photovoltaic and wind onshore, we can also implement offshore wind energy in a similar fashion.
Just use the profile you find for "Wind_onshore" in the profiles.csv file and multiply all its values with a factor of 1.2 and name the new column "Wind_offshore" (or whatever you fancy).

After you have added the new technology, check how the system costs have changed due to the addition. You can assume the following techno-economical data:

```
Life time: 20 years
Elec output per hour and unit: 1 GWh_el
Investment cost: 2550.0 million EUR/unit
Amortisation time: 25 years
Interest rate: 4.8 %
Fix operation and maintenance cost: 67.5 million EUR/year and unit
```

Also take a look at the dispatch behavior for the week with the highest / lowest total renewable energy generation (you can sum up the PV and Wind generation profiles before calculating the running average).

### (4) Analyse the effects of integration of renewable energies

Try to find out the contribution of the renewable technologies on the overall energy system. Look for the following indicators and fill the table:

- Renewable energy share (fraction of renewable generation to overall annual electricity demand)
- Carbon emission (in kt of CO2)
- Total system cost (in million EUR)

|                    | No renewable soures | PV only | Wind only | PV and Wind |
| ------------------ | ------------------- | ------- | --------- | ----------- |
| Renewable share    |                     |         |           |             |
| CO2 emissions      |                     |         |           |             |
| Total system costs |                     |         |           |             |

Hint: You can quickly disable individual technologies by setting the unitsUpperLimit to 0.
