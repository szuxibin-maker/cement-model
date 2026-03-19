from run.GEMSCalc import GEMS  # 你已有
import pandas as pd

def run_recipe_with_wc(oxide_row, verbose=True):
    from run.GEMSCalc import GEMS

    oxides_to_elements = {
        "SiO2": {"Si": 1, "O": 2},
        "Al2O3": {"Al": 2, "O": 3},
        "Fe2O3": {"Fe": 2, "O": 3},
        "CaO": {"Ca": 1, "O": 1},
        "MgO": {"Mg": 1, "O": 1},
        "K2O": {"K": 2, "O": 1},
        "Na2O": {"Na": 2, "O": 1},
        "SO3": {"S": 1, "O": 3},
        "CO2": {"C": 1, "O": 2},
    }

    oxide_molar_masses = {
        "SiO2": 60.08, "Al2O3": 101.96, "Fe2O3": 159.69, "CaO": 56.08,
        "MgO": 40.3, "SO3": 80.06, "CO2": 44.01,"K2O": 94.2, "Na2O": 61.98,
    }#"K2O": 94.2, "Na2O": 61.98,

    # 初始化 GEMS
    gem = GEMS("gems_files/MySystem2-dat.lst")
    gem.clear()
    # # gemsk.supress_species('CO2@')
    # gem.supress_species('zeoliteP_Ca')
    
    # # gem.supress_phase('ettringite-AlFe')
    # # gem.supress_phase('monosulph-AlFe')
    # gem.supress_phase('Al(OH)3mic')
    # gem.supress_phase('Gibbsite')
    # gem.supress_phase('Kaolinite')
    # gem.supress_phase('Graphite')
    # gem.supress_phase('CA')
    # gem.supress_phase('CA2')
    # # gem.supress_phase('C3AH6')
    # gem.supress_phase('C4AH13')
    # gem.supress_phase('C4AH19')
    # gem.supress_phase('CAH10')
    # # gem.supress_phase('C4AsH16')
    # gem.supress_phase('Chabazite')
    # gem.supress_phase('ZeoliteP')
    # # gem.supress_phase('Friedels')
    # # gem.supress_phase('Kuzels')
    # # gem.supress_phase('SO4_OH_AFm')
    # gem.supress_phase('Dolomite-dis')
    # gem.supress_phase('thaumasite')
    # gem.supress_phase('Hematite')
    # gem.supress_phase('Goethite')
    # gem.supress_phase('Pyrite')
    # gem.supress_phase('Magnesite')
    # gem.supress_phase('Pyrolusite')
    # gem.supress_phase('Natrolite')
    # gem.supress_phase('Quartz')
    # gem.supress_phase('ZeoliteX')
    # # gem.supress_phase('MSH')
    # gem.supress_phase('ZeoliteY')
    # # gem.supress_phase('ettringite')
    # # gem.supress_phase('Ferrihydrite-am')
    # # gem.supress_phase('ettringite-FeAl')
    # # gem.supress_phase('monosulph-FeAl')
    # # gem.supress_phase('OH_SO4_AFm')
    # # gem.supress_phase('CO3_SO4_AFt')
    # gem.supress_phase('hydrotalc-pyro')
    # gem.supress_phase('Al(OH)3am')
    # gem.supress_phase('Fe-carbonate')
    # # gem.supress_phase('Ferrihydrite-mc')
    # gem.supress_phase('Dolomite-ord')
    # gem.supress_phase('Rhodochrosite-sy')
    # gem.supress_phase('Ti(alpha)')
    # gem.supress_phase('TiO2(am_hyd)')
    # gem.supress_phase('Magnetite')
    # gem.supress_phase('syngenite')
     # gemsk.supress_species('CO2@')
    
    # gem.supress_species('zeoliteP_Ca')
    # # gem.supress_phase('ettringite-AlFe')
    # # gem.supress_phase('monosulph-AlFe')
    # # gem.supress_phase('Al(OH)3mic')
    # # gem.supress_phase('Gibbsite')
    # # gem.supress_phase('Kaolinite')
    # gem.supress_phase('Graphite')
    # gem.supress_phase('CA')
    # gem.supress_phase('CA2')
    # # gem.supress_phase('C3AH6')
    # # gem.supress_phase('C4AH13')
    # # gem.supress_phase('C4AH19')
    # # gem.supress_phase('CAH10')
    # # gem.supress_phase('C4AsH16')
    # gem.supress_phase('Chabazite')
    # gem.supress_phase('ZeoliteP')
    # # gem.supress_phase('Friedels')
    # # gem.supress_phase('Kuzels')
    # # gem.supress_phase('SO4_OH_AFm')
    # gem.supress_phase('Dolomite-dis')
    # # gem.supress_phase('thaumasite')
    # gem.supress_phase('Hematite')
    # gem.supress_phase('Goethite')
    # gem.supress_phase('Pyrite')
    # gem.supress_phase('Magnesite')
    # gem.supress_phase('Pyrolusite')
    # gem.supress_phase('Natrolite')
    # gem.supress_phase('Quartz')
    # gem.supress_phase('ZeoliteX')
    # # gem.supress_phase('MSH')
    # gem.supress_phase('ZeoliteY')
    # # gem.supress_phase('ettringite')
    # # gem.supress_phase('Ferrihydrite-am')
    # # gem.supress_phase('ettringite-FeAl')
    # # gem.supress_phase('monosulph-FeAl')
    # # gem.supress_phase('OH_SO4_AFm')
    # # gem.supress_phase('CO3_SO4_AFt')
    # gem.supress_phase('hydrotalc-pyro')
    # # gem.supress_phase('Al(OH)3am')
    # # gem.supress_phase('Fe-carbonate')
    # # gem.supress_phase('Ferrihydrite-mc')
    # gem.supress_phase('Dolomite-ord')
    # gem.supress_phase('Rhodochrosite-sy')
    # gem.supress_phase('Ti(alpha)')
    # gem.supress_phase('TiO2(am_hyd)')
    # # gem.supress_phase('Magnetite')
    # gem.supress_phase('syngenite')

        # === aq_gen ===
    # gem.supress_phase("aq_gen")
    
    # === gas_gen ===
    # gem.supress_phase("gas_gen")
    
    # === CSH ===
    # gem.supress_phase("CSHQ")
    
     # === AFt-phases（三钙铝石类，AFt） ===
    gem.supress_phase("ettringite-AlFe")     # 含铝铁的钙矾石，AFt 结构变体之一 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    gem.supress_phase("ettringite-FeAl")     # 含铁铝的钙矾石，另一种AFt结构混合型 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    # gem.supress_phase("ettringite")          # 钙矾石，C6AŜ3H32，典型的硫铝酸盐水化产物
    gem.supress_phase("SO4_CO3_AFt")         # 硫酸根-碳酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    gem.supress_phase("CO3_SO4_AFt")         # 碳酸根-硫酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    gem.supress_phase("C6AsH13")             # 硫酸型AFt矿物变体，可能水合度不同 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    gem.supress_phase("C6AsH9")              # 更低水合度的硫酸型AFt矿物 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    gem.supress_phase("thaumasite")          # 硅钙矾石（Thaumasite），Si 替代 SO₄ 形成的 AFt 型结构 需低温（<15°C）+高 CO₂ + Si源（如硅灰）+ SO₄²⁻，在 OPC 正常条件下极少生成，但在寒冷潮湿环境下可能导致腐蚀问题（称为“thaumasite 硫酸盐攻击”）
    
    # === AFm-Silica（硅铝酸盐型 AFm） ===
    # gem.supress_phase("straetlingite")       # 含硅 AFm，相当于 C2ASH8，常见于硅含量较高体系 在硅铝材料体系（如含火山灰、偏高岭土的胶凝材料）中可见，在纯 OPC 中较少见，但并非罕见
    
    # === AFm-phases（单盐型水化铝酸盐 AFm） ===
    gem.supress_phase("C4AH19")              # 高水合度铝酸钙水合物，早期水化产物 高水合度，早期生成，
    gem.supress_phase("C4AH13")              # 中等水合度 AFm，常在水化过程中生成 ，是常见稳定产物之一
    gem.supress_phase("C4AH11")              # 较低水合度的 AFm，相对更稳定 ，相对稳定，晚期可能形成
    gem.supress_phase("CAH10")               # 一种低温条件下形成的 AFm 类型 低温环境或高水灰比条件下更易生成，实验中难观察
    # gem.supress_phase("C4Ac0.5H12")          # 碳酸-氢氧共插层型 AFm，相对较高水合度（CO₃²⁻ + OH⁻）需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gem.supress_phase("C4Ac0.5H105")         # 同上，略低水合度版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gem.supress_phase("C4Ac0.5H9")           # 碳酸型 AFm 的低水合版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gem.supress_phase("C4AcH11")             # 纯碳酸插层型 AFm，水合度较高 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gem.supress_phase("C4AcH9")              # 纯碳酸插层型 AFm，水合度较低 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    gem.supress_phase("C4AsH16")             # 硫酸型 AFm，结构为 C4ASH16 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    gem.supress_phase("C4AsH14")             # 硫酸插层 AFm 的低水合版本 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gem.supress_phase("C4AsH12")             # C4ASH12，典型的单硫酸盐 AFm 相当于 monosulphate，常见于 OPC 与石膏反应后期
    gem.supress_phase("C4AsH105")            # 水合程度更低的硫酸插层型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    gem.supress_phase("C4AsH9")              # 极低水合度硫酸型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    gem.supress_phase("Friedels")            # Friedel’s 盐，Cl⁻ 插层型 AFm，相当于 C4AClH10 Cl⁻ 存在时极易生成，尤其在含氯体系或使用海砂水泥中
    gem.supress_phase("C2ASH55")             # 含硅 AFm 类型，C2ASH5.5，过渡相 低水合硅铝酸盐，硅含量较高或活性矿物体系中更常出现（如偏高岭土胶凝）
    gem.supress_phase("C4FH13")              # 含铁型 AFm，Fe 替代部分 Al Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    gem.supress_phase("C4Fc05H10")           # 混合铁铝型 AFm，水合度较低 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    gem.supress_phase("C4FcH12")             # Fe 和 Al 共掺，水合度较高 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    gem.supress_phase("monosulph-FeAl")      # 单硫酸盐型 AFm，含 Fe/Al 混合 Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    gem.supress_phase("monosulph-AlFe")      # 同类结构，不同比例 Fe/Al  Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    gem.supress_phase("SO4_OH_AFm")          # 硫酸根与 OH⁻ 共存插层型 AFm 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    gem.supress_phase("OH_SO4_AFm")          # OH⁻ 与硫酸根共存插层型 AFm，比例不同 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    gem.supress_phase("C2AH75")              # 低水合度的 C2AH 相，常为过渡态 中间水合物，水化过程中的过渡相，尤其在实验初期较常见
    gem.supress_phase("Kuzels")              # Kuzel’s 盐，Cl⁻ + SO₄²⁻ 双插层 AFm，特殊结构 Cl⁻+SO₄²⁻ 双插层型 AFm，较特殊，常见于海水或特殊氯盐体系

    
    # === Hydrogarnet（类水石榴石） ===
    gem.supress_phase("C3AH6")               # 水合铝酸钙石榴石，C₃AH₆，稳定的水化产物之一 是 OPC 水化中铝酸钙类相的重要稳定产物（六方 → 立方结构转变）有在 高温（> 60°C） 条件下，水合铝酸盐类才会从 AFt（如 ettringite）/AFm 转变为 C₃AH₆。
    gem.supress_phase("C3FH6")               # 铁代铝的类石榴石结构矿物，C₃FH₆ Fe 掺杂型类石榴石，需要高 Fe₂O₃，普通 OPC 中不易形成
    gem.supress_phase("C3FS0.84H4.32")       # 部分取代型含硅铁石榴石结构 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    gem.supress_phase("C3(AF)S0.84H")        # 铁铝硅共掺的水化产物 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    gem.supress_phase("C3FS1.34H3.32")       # 高硅含量的类石榴石矿物，结构更复杂 更复杂的高硅铁石榴石相，模拟中可能造成不必要相转移
    
    # === Al(OH)₃（铝氢氧化物） ===
    gem.supress_phase("Al(OH)3am")           # 非晶态氢氧化铝，模拟初期水化或胶凝态存在 初期水化中可能有非晶铝羟胶出现，但难以定量，压不压看模拟目标
    gem.supress_phase("Al(OH)3mic")          # 微晶氢氧化铝，具一定结晶度 微晶氢氧化铝，若系统中铝离子浓度较高，可析出；否则较少生成
    gem.supress_phase("Gibbsite")            # 三水铝石，Al(OH)₃，自然界稳定相，常见于铝土矿 常见于自然铝土矿环境中，不是 OPC 水化产物，除非有偏高岭土或极强碱性体系
    
    # === Hydrotalcite（类水滑石） ===
    gem.supress_phase("hydrotalc-pyro")      # 焙烧后再水化的水滑石（LDH）类矿物 焙烧型水滑石，需人工处理或Mg-Al来源，正常水泥中无来源
    gem.supress_phase("OH-hydrotalcite")     # 氢氧根型水滑石，层状结构，具阴离子吸附功能 镁硅水合物，模拟低Ca胶凝体系时才出现（如 MgO + SiO₂ 体系），OPC 中基本不生成
    
    # === M-S-H（镁硅水合物） ===
    gem.supress_phase("MSH")                 # 镁硅水合物，Mg-Si-H 凝胶态矿物，类似于低 Ca 的 C-S-H
    
    # === Zeolites（沸石类矿物） ===
    gem.supress_phase("ZeoliteP")            # 沸石 P，合成类沸石，常用于吸附/离子交换模拟 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gem.supress_phase("Natrolite")           # 钠沸石，Na₂Al₂Si₃O₁₀·2H₂O，自然沸石类型之一 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gem.supress_phase("Chabazite")           # 查巴茨沸石，具有较大孔径，阳离子交换能力强 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gem.supress_phase("ZeoliteX")            # 合成沸石 X，工业催化与分子筛常用 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gem.supress_phase("ZeoliteY")            # 沸石 Y，石化催化常用沸石，结构稳定 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。

    
    # === Clinkers（熟料矿物） ===
    # gem.supress_phase("Belite")         # 2CaO·SiO₂，低温稳定的硅酸二钙，相对水化慢
    # gem.supress_phase("Aluminate")      # 铝酸盐类，可能指 C₃A（三钙铝酸），对硫酸盐敏感
    # gem.supress_phase("Alite")          # 3CaO·SiO₂，高温形成，是水泥早期强度的主要来源
    # gem.supress_phase("Ferrite")        # 四钙铁铝酸盐，C₄AF，含铁熟料相
    gem.supress_phase("CA")             # 一钙铝酸盐，水化产物或熟料矿相 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gem.supress_phase("CA2")            # 二钙铝酸盐，水化速度较慢 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gem.supress_phase("lime")           # 氧化钙，游离石灰，可导致膨胀或不稳定  游离氧化钙（CaO），可能存在于部分高温水泥中，容易水化成 Portlandite。如果你不想模拟可能产生膨胀或碱反应，也可压制。
    gem.supress_phase("arcanite")       # K₂SO₄，硫酸钾，可能来自掺合料或污染源  这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gem.supress_phase("thenardite")     # Na₂SO₄，十水硫酸钠，高溶解度，可引起盐胀 这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gem.supress_phase("Na-oxide")       # 氧化钠，通常来自助熔剂或玻璃相成分 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    gem.supress_phase("K-oxide")        # 氧化钾，常出现在水泥熟料中，影响硅酸盐相稳定性 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    
    # === CaCO₃（碳酸钙类） ===
    gem.supress_phase("Aragonite")      # 霰石，CaCO₃ 的另一晶型，常在低温或有机环境中形成
    # gem.supress_phase("Calcite")        # 方解石，最稳定的 CaCO₃ 晶体，广泛存在于自然界 霰石为 CaCO₃ 的亚稳态晶型，需要低温或生物有机环境，OPC 水化中不易出现
    
    # === Portlandite（氢氧化钙） ===
    # gem.supress_phase("Portlandite")    # Ca(OH)₂，水泥水化的主要产物之一，提供碱性环境
    
    # === Sulfates（硫酸盐矿物） ===
    gem.supress_phase("Anhydrite")      # 无水硫酸钙，CaSO₄，常作为脱水后形式存在 若不使用脱水石膏作调凝剂可以压制，常见于特种水泥
    gem.supress_phase("Gypsum")         # 石膏，CaSO₄·2H₂O，常见调凝剂
    gem.supress_phase("hemihydrate")    # 半水石膏，CaSO₄·0.5H₂O，建筑石膏的主要形态 半水石膏主要是建筑石膏相，在 OPC 水泥中不是常见水化中间体

    
    # === Fe-phases（铁矿物） ===
    #原因：OPC 中 Fe₂O₃ 含量虽有，但远低于形成这些稳定铁矿相的需求。除非你使用钢渣、红泥、Fe 掺合料，否则这些相不会自然出现。
    gem.supress_phase("Iron")              # 金属铁，可能用于模拟原始还原性条件下的零价铁存在 
    gem.supress_phase("Fe-carbonate")      # 铁碳酸盐（可能为总称），例如菱铁矿等类矿物
    gem.supress_phase("Siderite")          # 菱铁矿，FeCO₃，典型的还原条件下形成的碳酸盐矿物
    gem.supress_phase("Hematite")          # 赤铁矿，Fe₂O₃，最稳定的三价铁氧化物，红褐色
    gem.supress_phase("Magnetite")         # 磁铁矿，Fe₃O₄，混合价态（Fe²⁺/Fe³⁺），具有强磁性
    gem.supress_phase("Ferrihydrite-am")   # 非晶铁水合氧化物，模拟新生成的纳米尺度氧化铁
    gem.supress_phase("Ferrihydrite-mc")   # 微晶型铁水合氧化物，Ferrihydrite 的结构化形式
    gem.supress_phase("Goethite")          # 针铁矿，α-FeOOH，常见于风化壤土和铁锈中
    gem.supress_phase("Pyrite")            # 黄铁矿，FeS₂，常见硫化物矿物，亮黄色金属光泽
    gem.supress_phase("Troilite")          # 斜方硫化亚铁，FeS，常见于陨石中，也存在于还原环境
    gem.supress_phase("Melanterite")       # 绿矾，FeSO₄·7H₂O，易溶于水，氧化还原环境中生成
    
    # === Mg-phases（镁矿物） ===
#OPC 中基本无 Mg 来源，除非有镁基掺合料（如 MgO、海砂、海水
    gem.supress_phase("Magnesite")         # 菱镁矿，MgCO₃，常见的镁碳酸盐矿物
    gem.supress_phase("Brucite")           # 氢氧化镁，Mg(OH)₂，在碱性环境中可析出，弱结晶性

    
    # === Mn-phases（锰矿物） ===
    gem.supress_phase("Rhodochrosite")       # 菱锰矿，MnCO₃，锰的主要碳酸盐矿物，粉红色晶体
    gem.supress_phase("Rhodochrosite-sy")    # 合成菱锰矿，模拟人工或非天然条件下形成的 MnCO₃
    gem.supress_phase("Hausmannite")         # 软锰矿，Mn₃O₄，混合价态氧化物，常见于风化环境中
    gem.supress_phase("Pyrolusite")          # 二氧化锰，MnO₂，天然最常见的锰氧化物，黑色矿物
    gem.supress_phase("Manganite")           # 水合氧化锰，MnO(OH)，一种含水氧化物，斜方晶系
    gem.supress_phase("Pyrochroite")         # 氢氧化锰，Mn(OH)₂，常见于还原环境或低氧条件下

    # === Other ===
    # === Clay minerals（粘土矿物） ===
    gem.supress_phase("Kaolinite")  # 高岭土，层状铝硅酸盐矿物，常见于土壤和沉积物中 常见于土壤或天然黏土系统，但在你的水泥+矿渣+粉煤灰+硅灰体系中不具形成条件。如果你没有添加天然黏土矿物，它不会生成
    
    # === Carbon phases（碳相） ===
    gem.supress_phase("Graphite")  # 石墨，碳的结晶形式，用于导电或润滑材料 石墨是纯碳的晶体，在任何正常水泥反应体系中都不会自然生成；除非模拟炭黑、还原环境，否则毫无意义
    
    # === Clinkers（熟料矿物） ===
    gem.supress_phase("Mayenite")  # 马叶矿，钙铝氧化物，在水泥熟料中形成于高温阶段 在高铝水泥或高温煅烧体系（如CSA或CA水泥）中可生成。你使用的是 OPC + 辅料体系，若未引入高铝熟料，则建议压制
    
    # === Carbonates（碳酸盐） ===
#白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gem.supress_phase("Dolomite-dis")  # 无序白云石，相较于有序白云石结构更松散 白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gem.supress_phase("Dolomite-ord")  # 有序白云石，CaMg(CO3)2，沉积岩中常见碳酸盐矿物 
    
    # === Sulfates（硫酸盐） ===
    gem.supress_phase("syngenite")  # 辛格石，K2Ca(SO4)2·H2O，一种钾钙复合硫酸盐矿物 K₂Ca(SO₄)₂·H₂O，特殊的钾钙硫酸盐，需高浓度 K⁺ + SO₄²⁻ 才可能生成，不属于正常 OPC 水化产物，可能带来模型不稳定或假阳离子迁移结果
    
    # === Elemental phases（元素相） ===
    gem.supress_phase("Sulphur")  # 硫，常以单质形式出现于火山或矿物环境中 单质硫只在极强还原或火山环境中形成，水泥体系无形成机制
    
    # === Silica phases（硅相） ===
    gem.supress_phase("Quartz")  # 石英，最常见的二氧化硅晶体结构 粉煤灰和矿渣中可能含有惰性石英，但通常不参与反应。若只关注活性反应产物，可压制；若要模拟总固相矿物平衡，可保留
    # gem.supress_phase("Silica-amorph")  # 非晶态二氧化硅，模拟水化产物或凝胶硅质成分 是硅灰、粉煤灰等活性成分形成的主要凝胶态结构之一，对 C-S-H（或 C-A-S-H）模拟至关重要，强烈建议保留
    
    # === Titanium oxides（钛氧化物） ===
#Ti 元素在你体系中基本不存在，TiO₂ 常用于颜料或光催化水泥，与结构矿物无关，压制无影响
    gem.supress_phase("Ti(alpha)")  # α-钛，六方结构的金属钛，在氧化之前的中间物相
    gem.supress_phase("TiO2(am_hyd)")  # 非晶水合二氧化钛，常见于纳米材料和胶体系统中

    


    
    # 计算总氧化物质量
    oxide_keys = list(oxides_to_elements.keys())
    oxide_mass_total = sum(oxide_row[oxide] for oxide in oxide_keys)

    # 初始水灰比
    wc_ratio = oxide_row["w/c"]
    water_mass = wc_ratio * oxide_mass_total  # 单位 g

    # 添加水
    gem.add_amt_from_formula({"H": 2, "O": 1}, water_mass * 1e-3, units="kg")

    # 添加氧化物
    for oxide in oxide_keys:
        mass = oxide_row[oxide]
        mol = mass / oxide_molar_masses[oxide]
        for el, ratio in oxides_to_elements[oxide].items():
            gem.add_element_amt(el, mol * ratio)

    gem.add_species_amt("O2@", 0.000001, units='kg')

    # === 尝试平衡 ===
    status = gem.equilibrate()
    if status.startswith("Failure") or status.startswith("Bad"):
        if verbose:
            print("⚠️ 初次反应失败，尝试加水再试一次...")
        
        # 再加 1g 水
        gem.add_amt_from_formula({"H": 2, "O": 1}, 0.001, units="kg")
        status_retry = gem.equilibrate()

        if status_retry.startswith("Failure") or status_retry.startswith("Bad"):
            if verbose:
                print("❌ 再次失败，跳过此配方")
            return None

        if verbose:
            print("✅ 加水后成功！")
            
    # === 提取结果 ===
    return {
        "b": gem.b.copy(),
        "pH": gem.pH,
        "phase_mass": gem.phase_masses,
        "phase_volume": gem.phase_volumes,
        "aq_composition": gem.aq_composition,
        "aq_element": gem.aq_element,
        "density": gem.density,
        "ionic_strength": gem.ionic_strength,
        "system_volume": gem.system_volume - gem.phase_volumes.get("gas_gen", 0.0),  # ✅ 修正这里
        "phase_sat_indices": gem.phase_sat_indices
    }
