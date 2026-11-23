# 序贯大坝代理模型的物理特征进化

你在为三个独立的坝体变形约束响应提出衍生物理特征。当前模型将物理特征的 Matérn-5/2 ARD 核与原始八维设计的 Matérn-5/2 ARD 残差核相加。物理候选用于几何表达，模型的拟合和采用由本地程序处理。

每项响应的 `matched_pairs` 包含五个绝对留一预测残差最大的断面，以及从低残差样本中为各断面匹配的对照。带符号残差为有限元真值减去留一预测，正值表示低估。响应和残差以米计。几何原语、原始变量及其别名由上下文定义。

1. 分别对比每项响应的高残差组与低残差对照，寻找可重复出现的几何区别。根据残差符号描述方向。符号混合时，可以在力学解释中说明适用的断面子组。
2. 每项响应给出恰好两条机制不同的新公式。每条解释应指出支持它的观察与可能的力学机制，注明这些解释是候选假设。
3. 公式只可使用上下文中的原始变量、别名、几何原语、数值常量和 `+`、`-`、`*`、`/`、`**`、`sqrt`、`abs`、`log`、`log1p`、`min`、`max`。幂指数为绝对值不大于 3 的常量。遵守量纲一致与函数定义域，避免分母在几何范围内触零。公式应组合多个变化量或对几何原语无量纲化。已有特征值用于诊断，表达式应直接写成原始变量与原语的组合。
4. 返回一个纯 JSON 对象，键为下面 `responses` 中的三个响应全名。每个键对应两元素数组，每个元素只有 `name`、`expression`、`physical_rationale` 三个字符串字段。`name` 使用英文 snake_case，且与该响应的已有特征及另一候选不同。公式只依赖设计几何，真值、预测、残差和样本编号用于诊断。

## 当前已观测前缀

```json
{
  "sample_count": 36,
  "response_unit": "m",
  "residual_definition": "observed minus LOO prediction",
  "original_variable_units": {
    "upstream_slope": "1",
    "downstream_slope": "1",
    "crest_width": "m",
    "modulus_zone_height": "m",
    "secondary_top_width": "m",
    "secondary_height": "m",
    "secondary_bottom_width": "m",
    "secondary_right_slope_distance": "m"
  },
  "aliases": {
    "x1": "upstream_slope",
    "x2": "downstream_slope",
    "x3": "crest_width",
    "x4": "modulus_zone_height",
    "x5": "secondary_top_width",
    "x6": "secondary_height",
    "x7": "secondary_bottom_width",
    "x8": "secondary_right_slope_distance",
    "d_ds_perp": "secondary_right_slope_distance"
  },
  "primitive_catalog": {
    "A_tot": {
      "unit": "m^2",
      "desc": "坝体主体域面积（GeoStudio 材料 4、5、6、12 的并集）"
    },
    "A_mod": {
      "unit": "m^2",
      "desc": "坝顶增模区(Material 12)面积（顶部高模量加强帽盖）"
    },
    "A_sec": {
      "unit": "m^2",
      "desc": "次堆石区(Material 5)面积（内部低模量软弱压缩区）"
    },
    "A_drain": {
      "unit": "m^2",
      "desc": "底部排水区(Material 6)面积（底部高刚度支承区）"
    },
    "A_main": {
      "unit": "m^2",
      "desc": "主堆石区(Material 4)面积（高刚度主承重与抗推骨架）"
    },
    "A_ds_shell": {
      "unit": "m^2",
      "desc": "次堆石区顶底高程之间、其右腰下游侧的材料 4 面积"
    },
    "B_bot": {
      "unit": "m",
      "desc": "坝体主体域在坝底高程的水平宽度"
    },
    "B_mod_bot": {
      "unit": "m",
      "desc": "坝顶增模区底面宽度"
    },
    "h_drain": {
      "unit": "m",
      "desc": "底部排水区高度（次堆石区底面至坝基的距离）"
    },
    "d_ds_perp": {
      "unit": "m",
      "desc": "次堆石区右腰到 GeoStudio 下游参考线的法向距离，等于 x8"
    },
    "d_us_min": {
      "unit": "m",
      "desc": "次堆石区左腰距主体域上游边界的最小水平净宽"
    },
    "X_c_sec": {
      "unit": "m",
      "desc": "次堆石区形心相对坝轴线的横坐标，正值偏下游"
    },
    "Y_c_sec": {
      "unit": "m",
      "desc": "次堆石区形心相对坝底的高程"
    },
    "W_sec_mean": {
      "unit": "m",
      "desc": "次堆石区平均水平宽度 0.5*(x5+x7)"
    },
    "H_dam": {
      "unit": "m",
      "desc": "大坝总高度常量"
    }
  },
  "responses": {
    "crest_max_settlement_m": {
      "current_formulas": [
        {
          "name": "narrow_crest_broad_base_index",
          "expression": "B_bot / crest_width",
          "physical_rationale": "候选假设：高残差断面相对低残差对照多数表现为坝顶更窄、坝底更宽，即更大的坝底宽度与坝顶宽度之比。正残差子组中这种几何可能使坝顶竖向荷载向更宽、更深的坝体扩散，放大坝顶下方压缩区而让模型低估沉降；负残差子组也出现该几何，说明方向会被帽盖和壳体刚度调节，单独使用时主要表征荷载扩散尺度。"
        },
        {
          "name": "soft_area_ratio",
          "expression": "A_sec / A_tot",
          "physical_rationale": "低模量次堆石占比表征主体断面的整体可压缩性。"
        },
        {
          "name": "elevated_drain_soft_cap_deficit",
          "expression": "(A_drain / A_tot) * (Y_c_sec / H_dam) * (A_sec / (A_mod + A_sec))",
          "physical_rationale": "候选假设：高绝对残差中的多数负残差断面以及一个正残差断面，相对对照表现为底部排水面积占比更大、软区形心更高、顶部帽盖相对软区更弱，因此该指标增大；其他正残差断面不完全支持。力学上，高位软区上方帽盖不足而底部排水刚体支承占比大，会改变竖向拱架与压缩柱路径：在排水支承占主导时模型可能把坝顶沉降算得过高，在帽盖极薄时又可能低估局部沉降，所以此解释仅作为符号混合子组的候选假设。"
        }
      ],
      "loo_nrmse": 0.22823387190384645,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 3,
        "overestimated": 2,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 7,
            "design_id": "lhs_007",
            "original_variables": {
              "upstream_slope": 1.3964213037264912,
              "downstream_slope": 1.55273578190388,
              "crest_width": 9.43626151664104,
              "modulus_zone_height": 12.050513005797814,
              "secondary_top_width": 12.720244476030498,
              "secondary_height": 20.88160543608512,
              "secondary_bottom_width": 35.83143862047909,
              "secondary_right_slope_distance": 13.627449547578786
            },
            "geometric_primitives": {
              "A_tot": 8042.333623728347,
              "A_mod": 294.3714524013453,
              "A_sec": 506.91854484567375,
              "A_drain": 3132.356073842461,
              "A_main": 4108.687552638867,
              "A_ds_shell": 522.9673958031354,
              "B_bot": 220.27388134400002,
              "B_mod_bot": 42.58935197700001,
              "h_drain": 38.06788155800001,
              "d_ds_perp": 13.627449547578786,
              "d_us_min": 4.945136433000016,
              "X_c_sec": 8.215963533000007,
              "Y_c_sec": 46.85203440499998,
              "W_sec_mean": 24.275841548499997,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 23.34334216527833,
              "soft_area_ratio": 0.06303127531914937,
              "elevated_drain_soft_cap_deficit": 0.1560034604301601
            },
            "loo_prediction_m": 0.08410383307917516,
            "observed_m": 0.07708631363929866,
            "observed_minus_loo_m": -0.007017519439876496
          },
          "matched_low_residual": {
            "sample_index": 32,
            "design_id": "eval_0032",
            "original_variables": {
              "upstream_slope": 1.3320104010403158,
              "downstream_slope": 1.3892031960189342,
              "crest_width": 7.83259486220777,
              "modulus_zone_height": 17.462718952447176,
              "secondary_top_width": 18.06828137487173,
              "secondary_height": 24.233660977333784,
              "secondary_bottom_width": 33.185473438352346,
              "secondary_right_slope_distance": 10.66262403037399
            },
            "geometric_primitives": {
              "A_tot": 7306.114505716314,
              "A_mod": 494.31824134870135,
              "A_sec": 621.0330589964988,
              "A_drain": 1779.5526228276067,
              "A_main": 4411.210582543508,
              "A_ds_shell": 438.5025185213277,
              "B_bot": 201.80239653600003,
              "B_mod_bot": 52.634710216,
              "h_drain": 29.303620070000022,
              "d_ds_perp": 10.66262403037399,
              "d_us_min": 16.62668582699999,
              "X_c_sec": 19.091129449999983,
              "Y_c_sec": 40.229172285000004,
              "W_sec_mean": 25.6268774065,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 25.76443695686286,
              "soft_area_ratio": 0.08500182395315617,
              "elevated_drain_soft_cap_deficit": 0.07372873442557262
            },
            "loo_prediction_m": 0.07286904797768441,
            "observed_m": 0.07437356009175546,
            "observed_minus_loo_m": 0.0015045121140710438
          }
        },
        {
          "high_residual": {
            "sample_index": 11,
            "design_id": "lhs_011",
            "original_variables": {
              "upstream_slope": 1.3717017078040479,
              "downstream_slope": 1.320912292949585,
              "crest_width": 9.763338568555401,
              "modulus_zone_height": 5.227712698589564,
              "secondary_top_width": 14.210190616724496,
              "secondary_height": 29.327369356524045,
              "secondary_bottom_width": 49.60905740491252,
              "secondary_right_slope_distance": 2.713831948671564
            },
            "geometric_primitives": {
              "A_tot": 7364.980022583283,
              "A_mod": 68.12974081013377,
              "A_sec": 935.825329374874,
              "A_drain": 2780.5265482474133,
              "A_main": 3580.498404150864,
              "A_ds_shell": 126.00460934764087,
              "B_bot": 201.616770115,
              "B_mod_bot": 20.240492020000005,
              "h_drain": 36.444917944999986,
              "d_ds_perp": 2.713831948671564,
              "d_us_min": 1.906027913999992,
              "X_c_sec": 16.229618721999998,
              "Y_c_sec": 48.39741523800001,
              "W_sec_mean": 31.909624011000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 20.650392148065347,
              "soft_area_ratio": 0.12706420472361732,
              "elevated_drain_soft_cap_deficit": 0.2301581931469669
            },
            "loo_prediction_m": 0.10135761338431643,
            "observed_m": 0.10697526562532986,
            "observed_minus_loo_m": 0.005617652241013432
          },
          "matched_low_residual": {
            "sample_index": 3,
            "design_id": "lhs_003",
            "original_variables": {
              "upstream_slope": 1.353045691022059,
              "downstream_slope": 1.508176082826362,
              "crest_width": 8.675030182722553,
              "modulus_zone_height": 6.283991266198675,
              "secondary_top_width": 17.933888328549898,
              "secondary_height": 31.60548164084404,
              "secondary_bottom_width": 69.30628530326103,
              "secondary_right_slope_distance": 2.4169107295218755
            },
            "geometric_primitives": {
              "A_tot": 7748.20338776067,
              "A_mod": 90.71643162123773,
              "A_sec": 1378.6338530328105,
              "A_drain": 3336.543297254967,
              "A_main": 2942.3098058516557,
              "A_ds_shell": 132.028256194689,
              "B_bot": 213.00543693799997,
              "B_mod_bot": 23.630378285999996,
              "h_drain": 33.11052709299997,
              "d_ds_perp": 2.4169107295218755,
              "d_us_min": 1.7045446040000058,
              "X_c_sec": 17.784611874000007,
              "Y_c_sec": 45.811391636999986,
              "W_sec_mean": 43.620086816500006,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 24.55385542775723,
              "soft_area_ratio": 0.17792948688086171,
              "elevated_drain_soft_cap_deficit": 0.2501273113506326
            },
            "loo_prediction_m": 0.10483866025666014,
            "observed_m": 0.10576399845433376,
            "observed_minus_loo_m": 0.0009253381976736252
          }
        },
        {
          "high_residual": {
            "sample_index": 19,
            "design_id": "lhs_019",
            "original_variables": {
              "upstream_slope": 1.6745232365833786,
              "downstream_slope": 1.4008319126769921,
              "crest_width": 9.165674889626958,
              "modulus_zone_height": 16.760990775170942,
              "secondary_top_width": 24.924136868277312,
              "secondary_height": 36.21510900519744,
              "secondary_bottom_width": 57.11954572428753,
              "secondary_right_slope_distance": 10.172838528352186
            },
            "geometric_primitives": {
              "A_tot": 8367.75443834113,
              "A_mod": 547.9646538156717,
              "A_sec": 1485.6104541145198,
              "A_drain": 1491.4065154909613,
              "A_main": 4842.7728149199775,
              "A_ds_shell": 626.5786803350254,
              "B_bot": 229.341951425,
              "B_mod_bot": 59.01079680599997,
              "h_drain": 18.023900219999973,
              "d_ds_perp": 10.172838528352186,
              "d_us_min": 16.988913295000003,
              "X_c_sec": 21.663753582000027,
              "Y_c_sec": 33.76287860899998,
              "W_sec_mean": 41.02184129599998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 25.02182918189172,
              "soft_area_ratio": 0.17753992006594252,
              "elevated_drain_soft_cap_deficit": 0.059407268784753164
            },
            "loo_prediction_m": 0.07921877804548015,
            "observed_m": 0.08383617259362591,
            "observed_minus_loo_m": 0.004617394548145762
          },
          "matched_low_residual": {
            "sample_index": 20,
            "design_id": "lhs_020",
            "original_variables": {
              "upstream_slope": 1.4420337151416434,
              "downstream_slope": 1.6952251044089008,
              "crest_width": 10.899606022680448,
              "modulus_zone_height": 13.154257709251336,
              "secondary_top_width": 23.49373414724033,
              "secondary_height": 32.683541574470254,
              "secondary_bottom_width": 52.227206782238454,
              "secondary_right_slope_distance": 7.27158587734832
            },
            "geometric_primitives": {
              "A_tot": 8659.800493040262,
              "A_mod": 386.1615117911777,
              "A_sec": 1237.4142604653896,
              "A_drain": 1956.4091925520615,
              "A_main": 5079.815528231634,
              "A_ds_shell": 461.85458236732734,
              "B_bot": 235.65675416000002,
              "B_mod_bot": 50.418226704000006,
              "h_drain": 25.162200715999973,
              "d_ds_perp": 7.27158587734832,
              "d_us_min": 12.972606106,
              "X_c_sec": 29.586810307000007,
              "Y_c_sec": 39.436926298,
              "W_sec_mean": 37.8604704645,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 21.62066717545878,
              "soft_area_ratio": 0.14289177463844333,
              "elevated_drain_soft_cap_deficit": 0.09176259820403383
            },
            "loo_prediction_m": 0.07629084062508133,
            "observed_m": 0.07761399443976472,
            "observed_minus_loo_m": 0.0013231538146833838
          }
        },
        {
          "high_residual": {
            "sample_index": 14,
            "design_id": "lhs_014",
            "original_variables": {
              "upstream_slope": 1.4751002650383611,
              "downstream_slope": 1.595669529967744,
              "crest_width": 8.064570860520767,
              "modulus_zone_height": 19.724240431447242,
              "secondary_top_width": 28.442154455937654,
              "secondary_height": 33.84352919365724,
              "secondary_bottom_width": 38.69247372526842,
              "secondary_right_slope_distance": 14.13000948719288
            },
            "geometric_primitives": {
              "A_tot": 8277.12235269257,
              "A_mod": 713.7342492409941,
              "A_sec": 1136.036374381157,
              "A_drain": 1186.1341314208164,
              "A_main": 5241.217597649602,
              "A_ds_shell": 893.9669002097403,
              "B_bot": 227.901531182,
              "B_mod_bot": 67.11122239899998,
              "h_drain": 17.432230375000017,
              "d_ds_perp": 14.13000948719288,
              "d_us_min": 12.44616542999998,
              "X_c_sec": 24.94246610099998,
              "Y_c_sec": 33.492772806000005,
              "W_sec_mean": 33.567314090500005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 28.25959807702444,
              "soft_area_ratio": 0.1372501608619572,
              "elevated_drain_soft_cap_deficit": 0.03983347277340295
            },
            "loo_prediction_m": 0.07526486563795061,
            "observed_m": 0.07099172116221743,
            "observed_minus_loo_m": -0.004273144475733182
          },
          "matched_low_residual": {
            "sample_index": 27,
            "design_id": "eval_0027",
            "original_variables": {
              "upstream_slope": 1.3056473664939403,
              "downstream_slope": 1.4340507615357638,
              "crest_width": 7.228919478133321,
              "modulus_zone_height": 14.876976828090847,
              "secondary_top_width": 10.501496735960245,
              "secondary_height": 30.423801876604557,
              "secondary_bottom_width": 45.867973673157394,
              "secondary_right_slope_distance": 6.909062474966049
            },
            "geometric_primitives": {
              "A_tot": 7313.767408192272,
              "A_mod": 361.42258919969555,
              "A_sec": 857.4867998220822,
              "A_drain": 1876.491797867288,
              "A_main": 4218.3662213032085,
              "A_ds_shell": 361.80544623935,
              "B_bot": 202.56657644299997,
              "B_mod_bot": 45.156990612000016,
              "h_drain": 25.699221294999973,
              "d_ds_perp": 6.909062474966049,
              "d_us_min": 22.954138701000005,
              "X_c_sec": 27.340317226000025,
              "Y_c_sec": 37.729782288000024,
              "W_sec_mean": 28.184735203999978,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 28.02169495119449,
              "soft_area_ratio": 0.11724283149360164,
              "elevated_drain_soft_cap_deficit": 0.09202675871878796
            },
            "loo_prediction_m": 0.0727176902090941,
            "observed_m": 0.07241404896544028,
            "observed_minus_loo_m": -0.0003036412436538266
          }
        },
        {
          "high_residual": {
            "sample_index": 15,
            "design_id": "lhs_015",
            "original_variables": {
              "upstream_slope": 1.330067888507769,
              "downstream_slope": 1.6538643380389764,
              "crest_width": 10.105264623039737,
              "modulus_zone_height": 8.088284809986758,
              "secondary_top_width": 16.701126741418904,
              "secondary_height": 21.205121572254274,
              "secondary_bottom_width": 55.28494324306234,
              "secondary_right_slope_distance": 6.607914647939461
            },
            "geometric_primitives": {
              "A_tot": 8184.754006807319,
              "A_mod": 156.67531762998547,
              "A_sec": 763.2366827567953,
              "A_drain": 4430.871938240758,
              "A_main": 2833.9700681797785,
              "A_ds_shell": 268.2284303453537,
              "B_bot": 223.51624487799998,
              "B_mod_bot": 31.70107872799997,
              "h_drain": 41.706593618,
              "d_ds_perp": 6.607914647939461,
              "d_us_min": 2.472503617000001,
              "X_c_sec": 11.089266749000018,
              "Y_c_sec": 50.41486541199998,
              "W_sec_mean": 35.9930349925,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 22.11879185908589,
              "soft_area_ratio": 0.09325102283123059,
              "elevated_drain_soft_cap_deficit": 0.3060013978084077
            },
            "loo_prediction_m": 0.08865321156696025,
            "observed_m": 0.09266916770421212,
            "observed_minus_loo_m": 0.004015956137251864
          },
          "matched_low_residual": {
            "sample_index": 2,
            "design_id": "lhs_002",
            "original_variables": {
              "upstream_slope": 1.569363259687807,
              "downstream_slope": 1.4665851746995955,
              "crest_width": 6.669684274972619,
              "modulus_zone_height": 7.144992267623804,
              "secondary_top_width": 15.054269429416003,
              "secondary_height": 22.4629175393885,
              "secondary_bottom_width": 51.30574787263373,
              "secondary_right_slope_distance": 5.155855220570123
            },
            "geometric_primitives": {
              "A_tot": 8083.082706026641,
              "A_mod": 106.02391874393064,
              "A_sec": 745.3197982798283,
              "A_drain": 3894.8369460715667,
              "A_main": 3336.9020429313186,
              "A_ds_shell": 202.66168602689874,
              "B_bot": 223.92986391,
              "B_mod_bot": 25.91713976700001,
              "h_drain": 41.392090193,
              "d_ds_perp": 5.155855220570123,
              "d_us_min": 1.969850661999999,
              "X_c_sec": 9.817879002000012,
              "Y_c_sec": 50.578356681,
              "W_sec_mean": 33.18000865150002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 33.57428248144764,
              "soft_area_ratio": 0.09220737005748161,
              "elevated_drain_soft_cap_deficit": 0.288325448359244
            },
            "loo_prediction_m": 0.09075505061704606,
            "observed_m": 0.09035779592611756,
            "observed_minus_loo_m": -0.00039725469092849464
          }
        }
      ]
    },
    "dam_max_settlement_m": {
      "current_formulas": [
        {
          "name": "soft_core_eccentricity",
          "expression": "X_c_sec / B_bot",
          "physical_rationale": "软区相对坝轴线的有符号偏心反映骨架承载分布的不对称性。"
        },
        {
          "name": "deep_soft_drainage_deficit",
          "expression": "((H_dam - Y_c_sec) / H_dam) * (A_sec / (A_sec + A_drain))",
          "physical_rationale": "候选假设：正残差断面的次堆石形心更低且 A_sec/(A_sec+A_drain) 更大，而多数负残差对照形心更高或底部排水支承更强；个别负残差断面因主骨架和帽盖很大而不服从该方向。力学上，低位软区若缺少底部高刚度排水支承，会形成更长的压缩柱并增大整体沉降，使模型低估。"
        }
      ],
      "loo_nrmse": 0.25269391202593927,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 2,
        "overestimated": 3,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 10,
            "design_id": "lhs_010",
            "original_variables": {
              "upstream_slope": 1.5468086093385431,
              "downstream_slope": 1.3012624657189515,
              "crest_width": 11.291715983318776,
              "modulus_zone_height": 13.534298401099962,
              "secondary_top_width": 26.23483907948985,
              "secondary_height": 27.391521289178613,
              "secondary_bottom_width": 63.350919984777626,
              "secondary_right_slope_distance": 8.592072960145197
            },
            "geometric_primitives": {
              "A_tot": 7897.915705416024,
              "A_mod": 372.63251196399375,
              "A_sec": 1226.9451132903607,
              "A_drain": 2947.1016370375714,
              "A_main": 3351.236443124097,
              "A_ds_shell": 381.73199744837643,
              "B_bot": 214.648971028,
              "B_mod_bot": 47.245817955999996,
              "h_drain": 30.074180309999974,
              "d_ds_perp": 8.592072960145197,
              "d_us_min": 7.236450769000015,
              "X_c_sec": 9.331825449999997,
              "Y_c_sec": 41.878520691000006,
              "W_sec_mean": 44.792879532000015,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.043474820332507914,
              "deep_soft_drainage_deficit": 0.12759440893394416
            },
            "loo_prediction_m": 0.1917529594215882,
            "observed_m": 0.17887868717138825,
            "observed_minus_loo_m": -0.012874272250199964
          },
          "matched_low_residual": {
            "sample_index": 1,
            "design_id": "lhs_001",
            "original_variables": {
              "upstream_slope": 1.6094092387020262,
              "downstream_slope": 1.398042551563892,
              "crest_width": 10.294596895348278,
              "modulus_zone_height": 17.22436933287436,
              "secondary_top_width": 22.75627311159097,
              "secondary_height": 28.75351770691302,
              "secondary_bottom_width": 46.6598200010875,
              "secondary_right_slope_distance": 12.621874746349105
            },
            "geometric_primitives": {
              "A_tot": 8262.411076914119,
              "A_mod": 581.5095874314916,
              "A_sec": 997.9784312431392,
              "A_drain": 1974.2684502131478,
              "A_main": 4708.654608026342,
              "A_ds_shell": 619.0967573010715,
              "B_bot": 225.44602486499997,
              "B_mod_bot": 60.22160053099998,
              "h_drain": 25.022112960000015,
              "d_ds_perp": 12.621874746349105,
              "d_us_min": 16.09770578199999,
              "X_c_sec": 15.802945496999996,
              "Y_c_sec": 37.74865155499998,
              "W_sec_mean": 34.70804655650002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.07009635901304095,
              "deep_soft_drainage_deficit": 0.16448591730430925
            },
            "loo_prediction_m": 0.16866897452355584,
            "observed_m": 0.1688936901258416,
            "observed_minus_loo_m": 0.00022471560228576326
          }
        },
        {
          "high_residual": {
            "sample_index": 19,
            "design_id": "lhs_019",
            "original_variables": {
              "upstream_slope": 1.6745232365833786,
              "downstream_slope": 1.4008319126769921,
              "crest_width": 9.165674889626958,
              "modulus_zone_height": 16.760990775170942,
              "secondary_top_width": 24.924136868277312,
              "secondary_height": 36.21510900519744,
              "secondary_bottom_width": 57.11954572428753,
              "secondary_right_slope_distance": 10.172838528352186
            },
            "geometric_primitives": {
              "A_tot": 8367.75443834113,
              "A_mod": 547.9646538156717,
              "A_sec": 1485.6104541145198,
              "A_drain": 1491.4065154909613,
              "A_main": 4842.7728149199775,
              "A_ds_shell": 626.5786803350254,
              "B_bot": 229.341951425,
              "B_mod_bot": 59.01079680599997,
              "h_drain": 18.023900219999973,
              "d_ds_perp": 10.172838528352186,
              "d_us_min": 16.988913295000003,
              "X_c_sec": 21.663753582000027,
              "Y_c_sec": 33.76287860899998,
              "W_sec_mean": 41.02184129599998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.09446049206171757,
              "deep_soft_drainage_deficit": 0.2713431235989366
            },
            "loo_prediction_m": 0.18585999955933752,
            "observed_m": 0.19813934916032805,
            "observed_minus_loo_m": 0.01227934960099053
          },
          "matched_low_residual": {
            "sample_index": 14,
            "design_id": "lhs_014",
            "original_variables": {
              "upstream_slope": 1.4751002650383611,
              "downstream_slope": 1.595669529967744,
              "crest_width": 8.064570860520767,
              "modulus_zone_height": 19.724240431447242,
              "secondary_top_width": 28.442154455937654,
              "secondary_height": 33.84352919365724,
              "secondary_bottom_width": 38.69247372526842,
              "secondary_right_slope_distance": 14.13000948719288
            },
            "geometric_primitives": {
              "A_tot": 8277.12235269257,
              "A_mod": 713.7342492409941,
              "A_sec": 1136.036374381157,
              "A_drain": 1186.1341314208164,
              "A_main": 5241.217597649602,
              "A_ds_shell": 893.9669002097403,
              "B_bot": 227.901531182,
              "B_mod_bot": 67.11122239899998,
              "h_drain": 17.432230375000017,
              "d_ds_perp": 14.13000948719288,
              "d_us_min": 12.44616542999998,
              "X_c_sec": 24.94246610099998,
              "Y_c_sec": 33.492772806000005,
              "W_sec_mean": 33.567314090500005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.10944404792559802,
              "deep_soft_drainage_deficit": 0.26779281991713155
            },
            "loo_prediction_m": 0.17248257667795708,
            "observed_m": 0.17191697423528282,
            "observed_minus_loo_m": -0.0005656024426742567
          }
        },
        {
          "high_residual": {
            "sample_index": 34,
            "design_id": "eval_0034",
            "original_variables": {
              "upstream_slope": 1.3704983707517386,
              "downstream_slope": 1.3009430952370167,
              "crest_width": 6.169699581339955,
              "modulus_zone_height": 11.456238799728453,
              "secondary_top_width": 15.937855262309313,
              "secondary_height": 20.547190830111504,
              "secondary_bottom_width": 25.968487756326795,
              "secondary_right_slope_distance": 3.5338027346879244
            },
            "geometric_primitives": {
              "A_tot": 7052.235868074875,
              "A_mod": 204.4022531363751,
              "A_sec": 430.52881350528736,
              "A_drain": 1618.647094475317,
              "A_main": 4798.657706957895,
              "A_ds_shell": 116.17872254426696,
              "B_bot": 196.45636355500002,
              "B_mod_bot": 33.516816694,
              "h_drain": 38.99657036999997,
              "d_ds_perp": 3.5338027346879244,
              "d_us_min": 12.045386390999994,
              "X_c_sec": 19.42913907000002,
              "Y_c_sec": 48.45047555399998,
              "W_sec_mean": 20.95317150999999,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.0988979879216823,
              "deep_soft_drainage_deficit": 0.07253942058683531
            },
            "loo_prediction_m": 0.14672879650593862,
            "observed_m": 0.15452222798622597,
            "observed_minus_loo_m": 0.007793431480287355
          },
          "matched_low_residual": {
            "sample_index": 11,
            "design_id": "lhs_011",
            "original_variables": {
              "upstream_slope": 1.3717017078040479,
              "downstream_slope": 1.320912292949585,
              "crest_width": 9.763338568555401,
              "modulus_zone_height": 5.227712698589564,
              "secondary_top_width": 14.210190616724496,
              "secondary_height": 29.327369356524045,
              "secondary_bottom_width": 49.60905740491252,
              "secondary_right_slope_distance": 2.713831948671564
            },
            "geometric_primitives": {
              "A_tot": 7364.980022583283,
              "A_mod": 68.12974081013377,
              "A_sec": 935.825329374874,
              "A_drain": 2780.5265482474133,
              "A_main": 3580.498404150864,
              "A_ds_shell": 126.00460934764087,
              "B_bot": 201.616770115,
              "B_mod_bot": 20.240492020000005,
              "h_drain": 36.444917944999986,
              "d_ds_perp": 2.713831948671564,
              "d_us_min": 1.906027913999992,
              "X_c_sec": 16.229618721999998,
              "Y_c_sec": 48.39741523800001,
              "W_sec_mean": 31.909624011000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.08049736494014263,
              "deep_soft_drainage_deficit": 0.08712244618321494
            },
            "loo_prediction_m": 0.1719364450609548,
            "observed_m": 0.1719201429519541,
            "observed_minus_loo_m": -1.6302109000693177e-05
          }
        },
        {
          "high_residual": {
            "sample_index": 23,
            "design_id": "eval_0023",
            "original_variables": {
              "upstream_slope": 1.3465167950838806,
              "downstream_slope": 1.3504548892378807,
              "crest_width": 9.308618433773518,
              "modulus_zone_height": 18.278893544338644,
              "secondary_top_width": 29.175809305161238,
              "secondary_height": 24.055996481329203,
              "secondary_bottom_width": 27.749729226343334,
              "secondary_right_slope_distance": 3.510688710026443
            },
            "geometric_primitives": {
              "A_tot": 7344.622640850989,
              "A_mod": 559.8062705784446,
              "A_sec": 684.7002773229959,
              "A_drain": 944.4711616301922,
              "A_main": 5155.644931319356,
              "A_ds_shell": 138.2280912578236,
              "B_bot": 201.48451856399998,
              "B_mod_bot": 55.86894898700001,
              "h_drain": 28.665109974000018,
              "d_ds_perp": 3.510688710026443,
              "d_us_min": 21.101340232000013,
              "X_c_sec": 29.066374193,
              "Y_c_sec": 40.793548714999986,
              "W_sec_mean": 28.4627692655,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.14426107971053515,
              "deep_soft_drainage_deficit": 0.18859251927106604
            },
            "loo_prediction_m": 0.15143471468806352,
            "observed_m": 0.1438072628529516,
            "observed_minus_loo_m": -0.007627451835111931
          },
          "matched_low_residual": {
            "sample_index": 17,
            "design_id": "lhs_017",
            "original_variables": {
              "upstream_slope": 1.6844315309886524,
              "downstream_slope": 1.4551043198867881,
              "crest_width": 6.360159107004971,
              "modulus_zone_height": 7.447526494976314,
              "secondary_top_width": 10.016459459486484,
              "secondary_height": 38.77411950955866,
              "secondary_bottom_width": 37.57497790629451,
              "secondary_right_slope_distance": 3.53078965948255
            },
            "geometric_primitives": {
              "A_tot": 8343.997175873314,
              "A_mod": 116.8886371461942,
              "A_sec": 922.65804003732,
              "A_drain": 1307.3043572125432,
              "A_main": 5997.146141477255,
              "A_ds_shell": 233.24576754904328,
              "B_bot": 231.285807562,
              "B_mod_bot": 27.62799936600001,
              "h_drain": 24.778353994999975,
              "d_ds_perp": 3.53078965948255,
              "d_us_min": 11.816037340999998,
              "X_c_sec": 32.28448626599999,
              "Y_c_sec": 40.423293504000014,
              "W_sec_mean": 23.795718683499985,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.13958697512101167,
              "deep_soft_drainage_deficit": 0.18773688403720884
            },
            "loo_prediction_m": 0.1432214915148393,
            "observed_m": 0.14316229337566488,
            "observed_minus_loo_m": -5.919813917440453e-05
          }
        },
        {
          "high_residual": {
            "sample_index": 35,
            "design_id": "eval_0035",
            "original_variables": {
              "upstream_slope": 1.6899565983563662,
              "downstream_slope": 1.464375388994813,
              "crest_width": 11.533354433253407,
              "modulus_zone_height": 15.195480976253748,
              "secondary_top_width": 27.14018240571022,
              "secondary_height": 37.962909899652004,
              "secondary_bottom_width": 31.129745049402118,
              "secondary_right_slope_distance": 2.0724369511008263
            },
            "geometric_primitives": {
              "A_tot": 8751.380236642724,
              "A_mod": 508.1272094412797,
              "A_sec": 1106.048002922161,
              "A_drain": 639.6401951867534,
              "A_main": 6497.564829092526,
              "A_ds_shell": 131.27636981507803,
              "B_bot": 237.553916987,
              "B_mod_bot": 57.89920367799999,
              "h_drain": 17.841609124,
              "d_ds_perp": 2.0724369511008263,
              "d_us_min": 27.507370690000016,
              "X_c_sec": 42.35913697699999,
              "Y_c_sec": 36.38986346799999,
              "W_sec_mean": 29.134963727500008,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.17831378035883141,
              "deep_soft_drainage_deficit": 0.3220182535949399
            },
            "loo_prediction_m": 0.15355054069023782,
            "observed_m": 0.14676447209809662,
            "observed_minus_loo_m": -0.006786068592141203
          },
          "matched_low_residual": {
            "sample_index": 4,
            "design_id": "lhs_004",
            "original_variables": {
              "upstream_slope": 1.5014740020205382,
              "downstream_slope": 1.5240172039007247,
              "crest_width": 8.143891063356346,
              "modulus_zone_height": 15.05385059674546,
              "secondary_top_width": 20.880637571229904,
              "secondary_height": 37.75439049568013,
              "secondary_bottom_width": 25.12606122906331,
              "secondary_right_slope_distance": 11.414546318035475
            },
            "geometric_primitives": {
              "A_tot": 8159.104536977323,
              "A_mod": 428.5172402523009,
              "A_sec": 868.4774359454105,
              "A_drain": 856.1050285227655,
              "A_main": 6006.004832256849,
              "A_ds_shell": 777.3720588287769,
              "B_bot": 224.63023579199998,
              "B_mod_bot": 51.727709488000016,
              "h_drain": 18.191758907999997,
              "d_ds_perp": 11.414546318035475,
              "d_us_min": 10.47389064699999,
              "X_c_sec": 28.690286130999993,
              "Y_c_sec": 36.488301820000004,
              "W_sec_mean": 23.003349400000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.12772228115170672,
              "deep_soft_drainage_deficit": 0.2552757608238603
            },
            "loo_prediction_m": 0.15979596540561813,
            "observed_m": 0.16155690587457133,
            "observed_minus_loo_m": 0.0017609404689531982
          }
        }
      ]
    },
    "dam_max_horizontal_displacement_m": {
      "current_formulas": [
        {
          "name": "soft_eccentric_area_fraction",
          "expression": "(X_c_sec / B_bot) * (A_sec / A_tot)",
          "physical_rationale": "将次堆石的相对面积与有符号偏心相乘以刻画软区对断面不对称的贡献。"
        },
        {
          "name": "main_drain_load_balance",
          "expression": "(A_main - A_drain) / A_tot",
          "physical_rationale": "候选假设：负残差高估子组中，高残断面相对对照表现为高刚度主堆石骨架面积相对底部排水支承更占优，即更大的主骨架与排水面积差占整体断面比例。主骨架占优可能把水平推力更多传给刚性承载路径，使实际水平位移小于模型预测；正残差和另一负残子组呈相反或较弱趋势，所以该解释仅作为特定刚度分配子组的候选假设。"
        },
        {
          "name": "slope_asymmetry_thrust_index",
          "expression": "(upstream_slope - downstream_slope) * (B_bot / H_dam)",
          "physical_rationale": "候选假设：高绝对残差断面多数具有更大的上游坡比相对下游坡比差，并与较大的坝底宽度耦合；正残差中的两个断面明显，负残差中的两个断面也如此，另一个正残差不服从。坡比不对称会移动自重合力与推力作用线，改变向下游剪切路径；当模型未充分捕捉这种几何偏心时，可能低估水平位移，而在主骨架或排水约束占优的负残差子组也可能转为高估，因此仅作为候选机制。"
        }
      ],
      "loo_nrmse": 0.24845136807814228,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 2,
        "overestimated": 3,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 19,
            "design_id": "lhs_019",
            "original_variables": {
              "upstream_slope": 1.6745232365833786,
              "downstream_slope": 1.4008319126769921,
              "crest_width": 9.165674889626958,
              "modulus_zone_height": 16.760990775170942,
              "secondary_top_width": 24.924136868277312,
              "secondary_height": 36.21510900519744,
              "secondary_bottom_width": 57.11954572428753,
              "secondary_right_slope_distance": 10.172838528352186
            },
            "geometric_primitives": {
              "A_tot": 8367.75443834113,
              "A_mod": 547.9646538156717,
              "A_sec": 1485.6104541145198,
              "A_drain": 1491.4065154909613,
              "A_main": 4842.7728149199775,
              "A_ds_shell": 626.5786803350254,
              "B_bot": 229.341951425,
              "B_mod_bot": 59.01079680599997,
              "h_drain": 18.023900219999973,
              "d_ds_perp": 10.172838528352186,
              "d_us_min": 16.988913295000003,
              "X_c_sec": 21.663753582000027,
              "Y_c_sec": 33.76287860899998,
              "W_sec_mean": 41.02184129599998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.016770508210026933,
              "main_drain_load_balance": 0.40050963781549614,
              "slope_asymmetry_thrust_index": 0.848228409632195
            },
            "loo_prediction_m": 0.10936074880633062,
            "observed_m": 0.12216219021832611,
            "observed_minus_loo_m": 0.01280144141199549
          },
          "matched_low_residual": {
            "sample_index": 9,
            "design_id": "lhs_009",
            "original_variables": {
              "upstream_slope": 1.521905918548933,
              "downstream_slope": 1.4341577509152494,
              "crest_width": 6.2300210381399586,
              "modulus_zone_height": 9.880374085179792,
              "secondary_top_width": 11.206422680062811,
              "secondary_height": 35.685761868782116,
              "secondary_bottom_width": 62.945727953476855,
              "secondary_right_slope_distance": 6.127392062091551
            },
            "geometric_primitives": {
              "A_tot": 7833.748071896116,
              "A_mod": 177.90726934915705,
              "A_sec": 1323.0879947750243,
              "A_drain": 2345.947228942543,
              "A_main": 3986.80557882939,
              "A_ds_shell": 374.8224797625352,
              "B_bot": 217.57872806899996,
              "B_mod_bot": 32.93087933900003,
              "h_drain": 25.433864045999997,
              "d_ds_perp": 6.127392062091551,
              "d_us_min": 11.431447804000015,
              "X_c_sec": 20.587154639999994,
              "Y_c_sec": 39.12681695700002,
              "W_sec_mean": 37.076075316499995,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.01598081857488956,
              "main_drain_load_balance": 0.2094601887662806,
              "slope_asymmetry_thrust_index": 0.2580018203259764
            },
            "loo_prediction_m": 0.08472689542399298,
            "observed_m": 0.08738836831204541,
            "observed_minus_loo_m": 0.002661472888052435
          }
        },
        {
          "high_residual": {
            "sample_index": 1,
            "design_id": "lhs_001",
            "original_variables": {
              "upstream_slope": 1.6094092387020262,
              "downstream_slope": 1.398042551563892,
              "crest_width": 10.294596895348278,
              "modulus_zone_height": 17.22436933287436,
              "secondary_top_width": 22.75627311159097,
              "secondary_height": 28.75351770691302,
              "secondary_bottom_width": 46.6598200010875,
              "secondary_right_slope_distance": 12.621874746349105
            },
            "geometric_primitives": {
              "A_tot": 8262.411076914119,
              "A_mod": 581.5095874314916,
              "A_sec": 997.9784312431392,
              "A_drain": 1974.2684502131478,
              "A_main": 4708.654608026342,
              "A_ds_shell": 619.0967573010715,
              "B_bot": 225.44602486499997,
              "B_mod_bot": 60.22160053099998,
              "h_drain": 25.022112960000015,
              "d_ds_perp": 12.621874746349105,
              "d_us_min": 16.09770578199999,
              "X_c_sec": 15.802945496999996,
              "Y_c_sec": 37.74865155499998,
              "W_sec_mean": 34.70804655650002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.008466615102115866,
              "main_drain_load_balance": 0.3309428848745256,
              "slope_asymmetry_thrust_index": 0.6439429649213041
            },
            "loo_prediction_m": 0.08989312003120725,
            "observed_m": 0.08021133624444485,
            "observed_minus_loo_m": -0.009681783786762405
          },
          "matched_low_residual": {
            "sample_index": 10,
            "design_id": "lhs_010",
            "original_variables": {
              "upstream_slope": 1.5468086093385431,
              "downstream_slope": 1.3012624657189515,
              "crest_width": 11.291715983318776,
              "modulus_zone_height": 13.534298401099962,
              "secondary_top_width": 26.23483907948985,
              "secondary_height": 27.391521289178613,
              "secondary_bottom_width": 63.350919984777626,
              "secondary_right_slope_distance": 8.592072960145197
            },
            "geometric_primitives": {
              "A_tot": 7897.915705416024,
              "A_mod": 372.63251196399375,
              "A_sec": 1226.9451132903607,
              "A_drain": 2947.1016370375714,
              "A_main": 3351.236443124097,
              "A_ds_shell": 381.73199744837643,
              "B_bot": 214.648971028,
              "B_mod_bot": 47.245817955999996,
              "h_drain": 30.074180309999974,
              "d_ds_perp": 8.592072960145197,
              "d_us_min": 7.236450769000015,
              "X_c_sec": 9.331825449999997,
              "Y_c_sec": 41.878520691000006,
              "W_sec_mean": 44.792879532000015,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.006753834853107899,
              "main_drain_load_balance": 0.0511698049409908,
              "slope_asymmetry_thrust_index": 0.7122463117275519
            },
            "loo_prediction_m": 0.09674769033282861,
            "observed_m": 0.09705632591519942,
            "observed_minus_loo_m": 0.0003086355823708059
          }
        },
        {
          "high_residual": {
            "sample_index": 34,
            "design_id": "eval_0034",
            "original_variables": {
              "upstream_slope": 1.3704983707517386,
              "downstream_slope": 1.3009430952370167,
              "crest_width": 6.169699581339955,
              "modulus_zone_height": 11.456238799728453,
              "secondary_top_width": 15.937855262309313,
              "secondary_height": 20.547190830111504,
              "secondary_bottom_width": 25.968487756326795,
              "secondary_right_slope_distance": 3.5338027346879244
            },
            "geometric_primitives": {
              "A_tot": 7052.235868074875,
              "A_mod": 204.4022531363751,
              "A_sec": 430.52881350528736,
              "A_drain": 1618.647094475317,
              "A_main": 4798.657706957895,
              "A_ds_shell": 116.17872254426696,
              "B_bot": 196.45636355500002,
              "B_mod_bot": 33.516816694,
              "h_drain": 38.99657036999997,
              "d_ds_perp": 3.5338027346879244,
              "d_us_min": 12.045386390999994,
              "X_c_sec": 19.42913907000002,
              "Y_c_sec": 48.45047555399998,
              "W_sec_mean": 20.95317150999999,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.006037579314488416,
              "main_drain_load_balance": 0.4509223276093657,
              "slope_asymmetry_thrust_index": 0.1846564391038974
            },
            "loo_prediction_m": 0.07965911880007723,
            "observed_m": 0.08870763001139202,
            "observed_minus_loo_m": 0.009048511211314794
          },
          "matched_low_residual": {
            "sample_index": 5,
            "design_id": "lhs_005",
            "original_variables": {
              "upstream_slope": 1.620861813812692,
              "downstream_slope": 1.6241883224247244,
              "crest_width": 7.396469826399812,
              "modulus_zone_height": 9.348096736641075,
              "secondary_top_width": 19.436622625497456,
              "secondary_height": 30.797597395220293,
              "secondary_bottom_width": 27.4739884022405,
              "secondary_right_slope_distance": 4.39242782746379
            },
            "geometric_primitives": {
              "A_tot": 8705.64496169005,
              "A_mod": 192.44218826179852,
              "A_sec": 722.3670559946428,
              "A_drain": 1236.183761906095,
              "A_main": 6554.651955527514,
              "A_ds_shell": 253.03518790104744,
              "B_bot": 240.130175399,
              "B_mod_bot": 36.05770447500001,
              "h_drain": 30.854305867999983,
              "d_ds_perp": 4.39242782746379,
              "d_us_min": 8.575461182999987,
              "X_c_sec": 29.596607445000018,
              "Y_c_sec": 45.37366048400003,
              "W_sec_mean": 23.455305514500026,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.010227092267275664,
              "main_drain_load_balance": 0.6109217889111952,
              "slope_asymmetry_thrust_index": -0.01079452833072475
            },
            "loo_prediction_m": 0.05702019585539838,
            "observed_m": 0.05880899645501255,
            "observed_minus_loo_m": 0.001788800599614171
          }
        },
        {
          "high_residual": {
            "sample_index": 31,
            "design_id": "eval_0031",
            "original_variables": {
              "upstream_slope": 1.3630816973745823,
              "downstream_slope": 1.3697291061282157,
              "crest_width": 6.585983991622925,
              "modulus_zone_height": 14.717251099646091,
              "secondary_top_width": 12.489079479128122,
              "secondary_height": 25.15388460829854,
              "secondary_bottom_width": 41.36039569508284,
              "secondary_right_slope_distance": 6.0032567055895925
            },
            "geometric_primitives": {
              "A_tot": 7249.350975337761,
              "A_mod": 343.7334144502847,
              "A_sec": 677.2617423696712,
              "A_drain": 2165.501912655367,
              "A_main": 4062.8539058624388,
              "A_ds_shell": 251.86419530972142,
              "B_bot": 201.41397894099998,
              "B_mod_bot": 43.944044393999974,
              "h_drain": 31.128864292000003,
              "d_ds_perp": 6.0032567055895925,
              "d_us_min": 21.600102883999995,
              "X_c_sec": 22.612904904999994,
              "Y_c_sec": 41.45810361100001,
              "W_sec_mean": 26.92473758699998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.010488760392196459,
              "main_drain_load_balance": 0.26172715318403666,
              "slope_asymmetry_thrust_index": -0.018092987117791184
            },
            "loo_prediction_m": 0.08559573920119998,
            "observed_m": 0.07656282041272364,
            "observed_minus_loo_m": -0.00903291878847634
          },
          "matched_low_residual": {
            "sample_index": 30,
            "design_id": "eval_0030",
            "original_variables": {
              "upstream_slope": 1.6313734091818333,
              "downstream_slope": 1.6791272994130848,
              "crest_width": 6.071929011493921,
              "modulus_zone_height": 19.592603002674878,
              "secondary_top_width": 14.385181814432144,
              "secondary_height": 20.009666811674833,
              "secondary_bottom_width": 61.474690763279796,
              "secondary_right_slope_distance": 6.464940689504147
            },
            "geometric_primitives": {
              "A_tot": 8790.413472151795,
              "A_mod": 725.9896508849653,
              "A_sec": 758.9653872971033,
              "A_drain": 3492.8760974358675,
              "A_main": 3812.5823365338565,
              "A_ds_shell": 250.80973158676687,
              "B_bot": 243.64897693799998,
              "B_mod_bot": 70.12193760200003,
              "h_drain": 31.39773018599999,
              "d_ds_perp": 6.464940689504147,
              "d_us_min": 43.306428987000004,
              "X_c_sec": 26.825576980999983,
              "Y_c_sec": 39.33241909899999,
              "W_sec_mean": 37.929936289000025,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.009505985480509556,
              "main_drain_load_balance": 0.036369874990616176,
              "slope_asymmetry_thrust_index": -0.15723224999532417
            },
            "loo_prediction_m": 0.053664925113820816,
            "observed_m": 0.05394199042999934,
            "observed_minus_loo_m": 0.0002770653161785208
          }
        },
        {
          "high_residual": {
            "sample_index": 2,
            "design_id": "lhs_002",
            "original_variables": {
              "upstream_slope": 1.569363259687807,
              "downstream_slope": 1.4665851746995955,
              "crest_width": 6.669684274972619,
              "modulus_zone_height": 7.144992267623804,
              "secondary_top_width": 15.054269429416003,
              "secondary_height": 22.4629175393885,
              "secondary_bottom_width": 51.30574787263373,
              "secondary_right_slope_distance": 5.155855220570123
            },
            "geometric_primitives": {
              "A_tot": 8083.082706026641,
              "A_mod": 106.02391874393064,
              "A_sec": 745.3197982798283,
              "A_drain": 3894.8369460715667,
              "A_main": 3336.9020429313186,
              "A_ds_shell": 202.66168602689874,
              "B_bot": 223.92986391,
              "B_mod_bot": 25.91713976700001,
              "h_drain": 41.392090193,
              "d_ds_perp": 5.155855220570123,
              "d_us_min": 1.969850661999999,
              "X_c_sec": 9.817879002000012,
              "Y_c_sec": 50.578356681,
              "W_sec_mean": 33.18000865150002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.004042697952430482,
              "main_drain_load_balance": -0.06902501476624247,
              "slope_asymmetry_thrust_index": 0.31101462951811676
            },
            "loo_prediction_m": 0.0673806140303409,
            "observed_m": 0.058805404918518496,
            "observed_minus_loo_m": -0.008575209111822406
          },
          "matched_low_residual": {
            "sample_index": 13,
            "design_id": "lhs_013",
            "original_variables": {
              "upstream_slope": 1.6446159216739957,
              "downstream_slope": 1.6795317296966799,
              "crest_width": 7.589823503377536,
              "modulus_zone_height": 11.310364914746476,
              "secondary_top_width": 27.938447008999383,
              "secondary_height": 24.297527356364657,
              "secondary_bottom_width": 65.86501425442128,
              "secondary_right_slope_distance": 8.027170857060547
            },
            "geometric_primitives": {
              "A_tot": 8935.356992263625,
              "A_mod": 279.49980368069726,
              "A_sec": 1139.5960830739239,
              "A_drain": 3871.8024894206746,
              "A_main": 3644.4586160883296,
              "A_ds_shell": 378.5524835054062,
              "B_bot": 246.17674519499997,
              "B_mod_bot": 43.878225472000025,
              "h_drain": 35.392107728999974,
              "d_ds_perp": 8.027170857060547,
              "d_us_min": 0.49380531299999575,
              "X_c_sec": 10.12346508600001,
              "Y_c_sec": 45.90354411200002,
              "W_sec_mean": 46.90173063150003,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.00524470666681819,
              "main_drain_load_balance": -0.025443177427514427,
              "slope_asymmetry_thrust_index": -0.1161548645253763
            },
            "loo_prediction_m": 0.04348652216592023,
            "observed_m": 0.04667452963675657,
            "observed_minus_loo_m": 0.003188007470836342
          }
        }
      ]
    }
  }
}
```
