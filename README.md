# Opening the Blackbox: Collision Attacks on Round-Reduced Tip5, Tip4, Tip4′ and Monolith

This repository accompanies the paper **Opening the Blackbox: Collision Attacks on Round-Reduced Tip5, Tip4, Tip4′ and Monolith**, accepted to ToSC 2024.4 (full version on ePrint). It contains the code used to produce the pratical collision attacks against round-reduced Tip5, Tip4, Tip4′ and Monolith.


## Repository structure

#### Targeted permutations:
- [`Tip5.sage`](Tip5.sage): This file contains a Sage implementation of the Tip5 permutation family, including [Tip5](https://eprint.iacr.org/2023/107) and [Tip4 & Tip4'](https://toposware.com/paper_tip5.pdf).
- [`Monolith.sage`](Monolith.sage): This file contains a Sage implementation of the [Monolith](https://tosc.iacr.org/index.php/ToSC/article/view/11810) permutation family.


#### Code to generate valid S-box differences:
- [`Diff_Property_of_Sbox_and_Finding_Diffs.cpp`](Diff_Property_of_Sbox_and_Finding_Diffs.cpp): This is the code to test the differential properties of the special S-boxes used in Tip5 and Monolith family. The code can also be used to search for the desired input-output differences for our efficient attacks.
- [`sys.sage`](sys.sage): This a script to generate the coefficient matrix 
(after Gaussian elimination) for the linear equation systems. 
The obtained matrices will be used in the file [`Diff_Property_of_Sbox_and_Finding_Diffs.cpp`](Diff_Property_of_Sbox_and_Finding_Diffs.cpp).


#### Algebraic attacks:
- [`collisionAttacks.ipynb`](collisionAttacks.ipynb): This file generates the equation systems for all practical attacks and creates the Magma scripts for running the Gröbner basis attack. 
- All Magma scripts are saved in the directory [`magma/`](magma/). If Gaussian elimination was performed on the equation system to reduce the number of equations and variables (as described in the paper), the attack scripts are saved in the [`magma/GE/`](magma/GE/) subdirectory. The output logs of running the collision attacks using Magma are saved in [`magma/results/`](magma/results/) and [`magma/GE/results/`](magma/GE/results/), respectively.
- [`results.ipynb`](results.ipynb): Summary and verification of the results.
- [`utils/`](utils/): Contains utility functions.




## Practical Results

The following practical results are sorted according to the corresponding chapters in the paper.


### 3.2 Collsion attack against 2-round Monolith-64
Magma files: [code](magma/GE/ica_mono64_t12r8c4d4.magma), [log](magma/GE/results/ica_mono64_t12r8c4d4.txt)


Note: This attack is a simple collision attack. It does not leverage the details of the S-Box.

| Part     |  Solution |
| -------- | --------- |
| IV       | `0000000000000000` `0000000000000000` `0000000000000000` `0000000000000000`          |
| $m_0$    | `0dafb1f39af4f126` `123c636bf95605db` `20d9d4075c10f47c` `9a1a6496a3fd807a` `44f906f66d366299` `db85d71b6755757f` `baba90fecb346771` `5080e621266de94d`          |
| $m_0'$   | `2855a6c059af53d7` `6c41bcca84872039` `548baaf10ba5cd31` `8ee31cb3244efa57` `a4df4c6bf38a6728` `0aaaacfc373fe34a` `902346c472ac3c2e` `5b865d4a1a3b69b4`          |
| Inner collision    | `000cda713a1183c5` `881fafb16f7a36a5` `7f3ce13953660b23` `d018d5dde86dcace`         |


### 5.3 Near collision attack against 3-round Tip4
Magma files: [code](magma/two_stage_nca_tip4_d2_s1_1_t1_4_u_1_m16r12c4d2_alpha7_pffffffff00000001.magma), [log](magma/results/two_stage_nca_tip4_d2_s1_1_t1_4_u_1_m16r12c4d2_alpha7_pffffffff00000001.txt)


| Part     | Solution |
| -------- | -------- |
| IV       | `0000000000000001` `0000000000000001` `0000000000000001` `0000000000000001`   |
| $m_0$    | `c6302a9416e1e7b3` `191740c2110b592d` `12b397fb06e33789` `29f060477c46d55f` `16994cb4ce4f0cc8` `076bdd4505177614` `44ffcea890f7f274` `295542230aec227e` `6cc6658831bc4f72` `2cc445c05ac6250c` `4721e7d54ca12911` `a487829e83cdc2a9`     |
| $m_0'$   | `1c0032b4ce793a22` `721d6b79e9cfb508` `c8b06bb2c12f1015` `a3458bf3e71280c6` `ab4623ee67400e0b` `4d3744b090c37448` `2425f22fe08cb9d4` `830f280903acb5ef` `9b62df4108d31f9c` `b63b522d993f6526` `3ed2820678c991ed` `06d101604089fa5c`    |
| Hash collision   | `ea4442b2c6cde046` `94587130d50f9581` `????????????????` `????????????????`  |




### 5.4 SFS collision attacks against 3-round Tip5, Tip4, Tip4'

#### Results for Tip5
Magma files: [code](magma/two_stage_sfs_ca_tip5_s1_0_t1_7_u_1_m16r10c6d5_alpha7_pffffffff00000001.magma), [log](magma/results/two_stage_sfs_ca_tip5_s1_0_t1_7_u_1_m16r10c6d5_alpha7_pffffffff00000001.txt)


| Part     | Solution |
| -------- | -------- |
| IV       | `63b5be0187fa05b0` `8f4db2bfcabd2302` `0a5c6021ee35d565` `11f40c2b7c8d452b` `bce318baebb56f57` `e98b25e991ed9829`         |
| $m_0$    | `0510053a63393b5e` `070b433153c77856` `a9df7ece7c71c3b6` `dea507128a64faa1` `a567a4cd0d87dbc1` `e4e0df5a55148b11` `189548e9e2d5578e` `f021f746d17c1717` `3c5e5bdf6e27f570` `f54e8c4543ed13c8`         |
| $m_0'$   | `ac95d195f2bd7273` `967bebc4a6922bf0` `4447e4bfc03b125f` `a8c7a57ed55a84ac` `76fda4ee213a7dc1` `de3f65a907005d3f` `ee35e7d8c8ae4ed0` `910d60b34fad5797` `c9821c0110240935` `201466ca5a0de071`         |
| Hash collision    | `acb46c393decc891` `e3abad309c016a4d` `2ffb1a0ba6264680` `4ba83933b09c17b5` `6fabfd5b7b5ac0e3`         |


#### Results for Tip4
Magma files: [code](magma/two_stage_sfs_ca_tip4_s1_0_t1_5_u_1_m16r12c4d4_alpha7_pffffffff00000001.magma), [log](magma/results/two_stage_sfs_ca_tip4_s1_0_t1_5_u_1_m16r12c4d4_alpha7_pffffffff00000001.txt)


| Part     | Solution |
| -------- | -------- |
| IV       | `05455d646a9fb498` `5b98a5cdcbf75625` `36611580308009e5` `f13961bc0bf79a35`         |
| $m_0$    |  `f065c68799b9224f` `5f6be20c98289f24` `2a35d4c999731105` `2f336da9b12d8e3b` `d5bf22405d1d1c6e` `62ea5eb72ef241af` `71dbff1985f169c4` `398e9ff03a7e62a8` `6ab41b1ef7261fee` `ae33d6bed1d6af17` `49f320563c9ad5fa` `91f85d88dd643c90`         |
| $m_0'$   | `f99251c7d59303fc` `0244d3ac998244b2` `5071809bfcedd711` `d2587fef54685ba9` `bf83eeefc8aafcd7` `b2a6864b9ead20c0` `e6299df48ae58aba` `1fb4938ab0faf420` `53d74eab02fe04a9` `f690f3967d1c28c8` `e0e56a7b656bda5d` `8a57287029e36c73`         |
| Hash collision    | `199feda6f6577448` `73977eac80d5f0db` `7864fb871040f78a` `08ac2e5ec2d03347`         |


#### Results for Tip4'
Magma files: [code](magma/two_stage_sfs_ca_tip4p_s1_0_t1_5_u_1_m12r8c4d4_alpha7_pffffffff00000001.magma), [log](magma/results/two_stage_sfs_ca_tip4p_s1_0_t1_5_u_1_m12r8c4d4_alpha7_pffffffff00000001.txt)

| Part     | Solution |
| -------- | -------- |
| IV       | `84aefeb85e1c736a` `99852f9ac85c8be0` `aa2ccabaf5e705e1` `6206827be78f9cd2`         |
| $m_0$    | `f9e8311f95cb2b77` `4835b953ea03ed41` `0b47521c88b06001` `ac8436d8c449aa10` `4ae194055c249e8e` `41ffa5ff7fc583ca` `dd37c9a863fe195b` `ad014cca15abdc41`         |
| $m_0'$   | `abf3f56fbab88503` `2f73090e75f87967` `1144a30cef40aee0` `c975f46f94e8a65e` `22b84ee8dd304de5` `b37f710c695a5ea3` `bc57a1d2c071fd0a` `16604e2c176495ae`         |
| Hash collision    | `4f885efd35b10b46` `3ce80cda16752fc8` `93a563d23acb5028` `a2cc2fdf2b4b0d6f`         |


### 7.2 Collsion attack against 2-round Monolith-31
Magma files: [code](magma/two_stage_ica_mono31_t24r16c8d8.magma), [log](magma/results/two_stage_ica_mono31_t24r16c8d8.txt)

| Part     |  Solution |
| -------- | --------- |
| IV       | `00000000` `00000000` `00000000` `00000000` `00000000` `00000000` `00000000` `00000000`          |
| $m_0$    | `29c00f10` `3b7173a8` `79ea0cd6` `0d55fe73` `04b4aaef` `578ce0fd` `1945b9f5` `5178fbea` `6d2b218c` `7942abf3` `583952bf` `6185fe9f` `5339a47f` `233560f3` `0010614c` `5b831c47`          |
| $m_0'$   | `622558a9` `3751266d` `7586a6a1` `6190bb51` `43728d5a` `43d555ae` `299fa6fa` `55327b56` `6ba000e3` `12f56ade` `4c930fde` `5b16a551` `5bd241bf` `4fc1844e` `47adf609` `2c7a8474`         |
| Inner collision    | `7c957753` `09b10941` `0bdfbfee` `51f1a150` `56e2ccfb` `31650e10` `5855422f` `67e6d6af`         |

