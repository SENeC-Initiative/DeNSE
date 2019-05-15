Remarks for cleanup and merge requests
======================================


Structure
---------

### Restructuring

* If we keep the `netgrowth` folder, move `Readme`, `docs` and `cmake` outside,
  rename `netgrowth` as`src`.
* create `data_io.py` that would regroup `data_to_json` and `swc_detection`?
* create testing folder and move `TestRandomGen` there


Code
----

**`CreateNeurons`**

   * Explicit is better than implicit, restore keywords arguments for
     `CreateNeurons`, remove `culture`: it is set through `CreateEnvironment`,
     here the objects are just interfering.
   * To simply things we can probably also save the corresponding culture
     inside the Python package.
   * keep the `culture` keyword but check it's contained in the main culture

**`Branching`**

	* l77 empty if
	* update me on deprecated stuff: `branching_event` vs `initialize_next_event` or `compute_next_event`


Models
------

**Lurd:**

* utilized to used
* initialize left = received
* initialize_T_distributions needs to be renamed (initialize_cr_distribution?)
* idem: `compute_cr_demand` (replace cr everywhere)
* compute_critical_resource into update_cr?
* what is `initialize_T_distributions` for?
* l116: ``critical_.demand= 1+0.1*prev_demand_;`` why the magic number?
* update the docstring (not multiplicative anymore)
* what's the plan for `after_split`?
* l157: why the floor?
* l178: don't use auto + add comment remind people what Branching is
* l185: precompute the pow and update it when necessary?
* `get_critical_resource_quotient` is unclear
* l193-196: not sure what's happening here...


Testing
-------

* After correction to make it work, ``__growth_test__.py`` burns up my RAM in
  5 seconds for 10 neurons. The problem starts with 6 neurons
* Plot make neurite start very far away from soma.
* Neuron make 180° turns + always at the same position
* `f_coeff` is wrong
* `move_` may be multiply defined between child classes and GC
* check that seeds is list
* does not seem to work with btmorph when several neurons are present


Documentation
-------------

* Links are with backquotes and not `'`
* Inverse the order of the recursive algorithm.
