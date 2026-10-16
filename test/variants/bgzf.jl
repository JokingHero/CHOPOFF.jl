using CodecZlib, Random

@testset "VCF BGZF reader" begin
    mktempdir() do dir
        rng = MersenneTwister(7)
        # Several batches of blocks, lines longer than one block, and CRLF ends.
        lines = [randstring(rng, "ACGT\t|/:", rand(rng, [0, 1, 50, 3000, 70000])) for _ in 1:400]
        text = join((l * (isodd(i) ? "\n" : "\r\n") for (i, l) in enumerate(lines)))
        bgzf = joinpath(dir, "a.vcf.gz")
        stream = CHOPOFF.VariantCallFormat.BCF.BGZFStreams.BGZFStream(bgzf, "w")
        write(stream, text)
        close(stream)
        gzip = joinpath(dir, "b.vcf.gz")
        write(gzip, transcode(GzipCompressor, text))
        read_all(path; keep = false) = (s = CHOPOFF.open_variant_vcf(path); try collect(eachline(s; keep)) finally close(s) end)
        @test CHOPOFF.open_variant_vcf(bgzf) isa CHOPOFF.VariantBGZFStream
        @test read_all(bgzf) == lines
        @test join(read_all(bgzf; keep = true)) == text
        @test read_all(gzip) == lines

        bytes = read(bgzf)
        corrupt = joinpath(dir, "corrupt.vcf.gz")
        second = CHOPOFF.variant_bgzf_block_size(bytes[1:18]) + 1
        bytes[second + 40] ⊻= 0xff
        write(corrupt, bytes)
        @test_throws TaskFailedException read_all(corrupt)
        truncated = joinpath(dir, "truncated.vcf.gz")
        write(truncated, read(bgzf)[1:(second + 100)])
        @test_throws TaskFailedException read_all(truncated)
    end
end
