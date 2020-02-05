using DataFrames, CSV, Gadfly, Compose, TableView
filterminlog2foldchange = 0.50
filterminTPM = 2.5

dst = CSV.read("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/extractedData/DeSeqOutputAllConds.tsv")
metadata = CSV.read("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/extractedData/sampleMetadata_SI2-SI4.txt")

######## create a series of "and" filters from singleton conditions ############
dstboolhigh = dst[:, r".*high_isDeGene"]
boolselectionvec1 = dstboolhigh[:, Symbol("RA-high_isDeGene")] .&
                    dstboolhigh[:, Symbol("TGFb-high_isDeGene")]

dstlogfoldchangehigh = dst[:, r".*high_log2fc"]
boolselectionvec2 = (dstlogfoldchangehigh[:, Symbol("RA-high_log2fc")] .> filterminlog2foldchange) .&
                    (dstlogfoldchangehigh[:, Symbol("TGFb-high_log2fc")] .> filterminlog2foldchange)

dstavgtpmhigh = dst[:, r".*high_avgTPM"]
boolselectionvec3 = (dstavgtpmhigh[:, Symbol("RA-high_avgTPM")] .> filterminTPM) .&
                    (dstavgtpmhigh[:, Symbol("TGFb-high_avgTPM")] .> filterminTPM)

finalselectionvec = boolselectionvec1 .& boolselectionvec2 .& boolselectionvec3
println(sum(finalselectionvec))
finalselectionbool = convert(Array{Bool, 1}, finalselectionvec)
################################################################################

dstfiltered = dst[finalselectionbool, :]
dstfiltselect = dstfiltered[:, r"hgnc.symbol|.*_log2fc"]
dststacked = stack(dstfiltered, r".*_log2fc", variable_name=:condition, value_name=:log2fc)
dststacked[:, :foldchange] = 2 .^ dststacked[:, :log2fc]

# now build long table with condition column
dst_tpms = dst[:, r".*_tpm|ensg|hgnc.*"]
dst_tpms_filt = dst_tpms[finalselectionbool, :]
dst_tpms_stacked = stack(dst_tpms_filt, r".*tpm", variable_name=:sampleID, value_name=:tpm)
conditions = [split(join(split(string(s), "-")[2:end], "-"), "_")[1] for s in dst_tpms_stacked.sampleID]
dst_tpms_stacked[:, :condition] = conditions

# add replicate info
dst_tpms_stacked[:, :sampleID] = map(s -> split(string(s), "_")[1], dst_tpms_stacked.sampleID)
dst_tpms_stacked = join(dst_tpms_stacked, metadata[:, [:sampleID; :replicate]], on = :sampleID, kind = :left)

function lessthanfx(left1, right1)
    left = string(left1)
    right = string(right1)
    if (left[1:3] != right[1:3]) | (occursin("and", left) ⊻ occursin("and", right))
        if (occursin("and", left))
            return false
        elseif (occursin("and", right))
            return true
        else
            return isless(left1, right1)
        end
    elseif occursin("low", left)
        return true
    elseif occursin("high", left)
        return false
    elseif occursin("med", left)
        return occursin("high", right)
    end
    return false
end

sort!(dststacked, :condition, lt=lessthanfx)
sort!(dst_tpms_stacked, :condition, lt=lessthanfx)
dstfiltered_avgTPM = stack(dstfiltered[:, r".*_avgTPM|ensg"], r".*_avgTPM", variable_name = :condition_tpm, value_name = :avg_tpm)
sort!(dstfiltered_avgTPM, :condition_tpm, lt=lessthanfx)

a = []
b = []
c = []
for (i,j) in enumerate(unique(dststacked.ensg))
    tempfiltereddst = dst_tpms_stacked[dst_tpms_stacked.ensg .== j, :]
    hgnc = tempfiltereddst[1, Symbol("hgnc.symbol")]
    println(hgnc)

    tempFilteredDstAvgTpm = dstfiltered_avgTPM[dstfiltered_avgTPM.ensg .== j, :]

    # calculate average variance of the "both" conditions
    variances = []
    for cond in ["TGFb-and-RA-low"; "TGFb-and-RA-med"; "TGFb-and-RA-high"]
        tpms = tempfiltereddst[tempfiltereddst.condition .== cond, :tpm]

    end

    # calculate predictions for additive and multiplicative models
    control_tpm_value = tempFilteredDstAvgTpm[tempFilteredDstAvgTpm.condition_tpm .== Symbol("EtOH-nlDensity_avgTPM"), Symbol("avg_tpm")]
    control_tpm_value = control_tpm_value[1]
    prediction_df = DataFrame(ensg = String[],
                              dosage = String[],
                              additive_pred_TPM = Float64[],
                              multiplicative_pred_TPM = Float64[])
    for dosage in ["low"; "med"; "high"]
        ra_tpm   = tempFilteredDstAvgTpm[tempFilteredDstAvgTpm.condition_tpm .== Symbol("RA-"*dosage*"_avgTPM"), Symbol("avg_tpm")][1]
        tgfb_tpm = tempFilteredDstAvgTpm[tempFilteredDstAvgTpm.condition_tpm .== Symbol("TGFb-"*dosage*"_avgTPM"), Symbol("avg_tpm")][1]
        ra_effect = ra_tpm - control_tpm_value
        tgfb_effect = tgfb_tpm - control_tpm_value
        addpred  = control_tpm_value + ra_effect + tgfb_effect
        multpred = (ra_tpm/control_tpm_value) * (tgfb_tpm/control_tpm_value) * control_tpm_value
        push!(prediction_df, (j, dosage, addpred, multpred))

        #for beeswarm plot without layers working
        push!(tempfiltereddst, ("NA", addpred, j, hgnc, "TGFb-and-RA-"*dosage, "addpred"))
        push!(tempfiltereddst, ("NA", multpred, j, hgnc, "TGFb-and-RA-"*dosage, "multpred"))
    end

    print(prediction_df)
    ## attempt to overlay lines on beeswarm plot seems to keep ignoring "x axis" values...
    beeswarmplot = layer(x = tempfiltereddst.condition,
                         y = tempfiltereddst.tpm,
                         color = tempfiltereddst.replicate,
                         Geom.beeswarm)
    lineplot_add_pred_low = layer(x = [1;2;3;0], y = [prediction_df[:, 3]; prediction_df[3, 3]], Geom.line)
    p = plot(beeswarmplot, lineplot_add_pred_low, Guide.title(hgnc))
    if hgnc == "NA"
        hgnc = j
    end

    push!(b, tempfiltereddst)
    push!(c, tempFilteredDstAvgTpm)

    # p = plot(tempfiltereddst, x = :condition, y = :tpm, color = :replicate, Geom.beeswarm, Guide.title(hgnc))
    draw(SVG(joinpath("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/plots/tpmBeeSwarmsAll", hgnc * ".svg")), p)
    push!(a, p)
    break
end
a[1]
# showtable(b[1])


# plot(dststacked, x="condition", y="log2fc", Geom.beeswarm)
# gridstack([a[1] a[2]; a[3] a[4]])
